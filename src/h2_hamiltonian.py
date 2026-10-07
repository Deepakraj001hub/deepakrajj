"""
H2 Molecular Hamiltonian Generator
====================================
Qiskit Fall Fest 2026 – Challenge I5: Molecule Energy Estimator

SOURCE & CONVENTIONS
---------------------
Coefficients are taken from Table I of:

    O'Malley, P. J. J., et al.
    "Scalable Quantum Simulation of Molecular Energies."
    Physical Review X 6, 031007 (2016).
    arXiv: 1512.06860v2
    DOI: 10.1103/PhysRevX.6.031007

    NOTE: The v2 arXiv revision corrects an error in the identity
    coefficient (g_I) of the original v1 Table I.
    This module uses the corrected v2 values.

HAMILTONIAN FORM
-----------------
The 2-qubit Bravyi-Kitaev (BK) reduced Hamiltonian for molecular hydrogen
in the minimal basis of Hartree-Fock orbitals, represented by the
O'Malley et al. Table I coefficients:

    H(R) = g_I(R) * II
           + g_Z0(R) * ZI
           + g_Z1(R) * IZ
           + g_Z0Z1(R) * ZZ
           + g_Y0Y1(R) * YY
           + g_X0X1(R) * XX

Qubit register convention (Qiskit little-endian, qubit 0 is rightmost):
    - 'II'  = I ⊗ I
    - 'ZI'  = Z ⊗ I   (Z on qubit 1, I on qubit 0)
    - 'IZ'  = I ⊗ Z   (I on qubit 1, Z on qubit 0)
    - 'ZZ'  = Z ⊗ Z
    - 'YY'  = Y ⊗ Y
    - 'XX'  = X ⊗ X

BOND DISTANCE UNITS
--------------------
All distances in this module are in Ångströms (Å).
The original Table I reports distances in Ångströms.

ENERGY UNITS
-------------
All energies are in Hartree (atomic units).
The identity coefficient g_I includes the nuclear repulsion energy
1/R (in atomic units, with R in Bohr; the conversion is performed
internally) plus constant electronic terms from the BK mapping.

REPRODUCIBILITY
----------------
Published tabulated coefficients are maintained in h2_table_i_data.py.
Cubic spline interpolation (scipy.interpolate.CubicSpline, natural
boundary conditions) is used only between the published table rows,
within the supported range [0.20, 2.85] Å. Interpolated values are not
published table entries.

SUPPORTED RANGE
---------------
R_MIN = 0.20 Å  (first published table point)
R_MAX = 2.85 Å  (last published table point)

Only requests within [R_MIN, R_MAX] are served; values outside this
range raise ValueError.

DEPENDENCIES
------------
    qiskit >= 2.5.2          (SparsePauliOp)
    numpy  >= 2.0            (array operations)
    scipy  >= 1.15           (CubicSpline, linalg.eigh)
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh
from qiskit.quantum_info import SparsePauliOp

from h2_table_i_data import COEFFICIENT_NAMES, PAULI_ORDER, TABLE_I

# Convert source columns (R, identity, Z0, Z1, Z0Z1, X0X1, Y0Y1)
# explicitly to internal order (R, g_I, g_Z0, g_Z1, g_Z0Z1, g_Y0Y1, g_X0X1).
_TABLE_I = np.array([
    (row[0], row[1], row[2], row[3], row[4], row[6], row[5])
    for row in TABLE_I
])

# Supported interpolation range
R_MIN: float = float(_TABLE_I[0, 0])   # 0.20 Å
R_MAX: float = float(_TABLE_I[-1, 0])  # 2.85 Å

# Build cubic splines once at module load time (natural BC)
_R_POINTS = _TABLE_I[:, 0]
_SPLINES = {
    name: CubicSpline(_R_POINTS, _TABLE_I[:, index + 1], bc_type="natural")
    for index, name in enumerate(COEFFICIENT_NAMES)
}


def _validate_distance(bond_distance_angstrom: float) -> None:
    """Raise ValueError if bond_distance_angstrom is outside [R_MIN, R_MAX]."""
    if not (R_MIN <= bond_distance_angstrom <= R_MAX):
        raise ValueError(
            f"Bond distance {bond_distance_angstrom:.4f} Å is outside the "
            f"supported interpolation range [{R_MIN:.2f}, {R_MAX:.2f}] Å "
            f"(O'Malley et al., PRX 6, 031007, 2016, Table I)."
        )


def interpolated_coefficients(bond_distance_angstrom: float) -> dict:
    """
    Return the six interpolated BK Hamiltonian coefficients at a given
    bond distance.

    Parameters
    ----------
    bond_distance_angstrom : float
        H–H internuclear distance in Ångströms.

    Returns
    -------
    dict with keys:
        'g_I', 'g_Z0', 'g_Z1', 'g_Z0Z1', 'g_Y0Y1', 'g_X0X1'
    All values are floats in Hartree.

    Raises
    ------
    ValueError
        If bond_distance_angstrom is outside [R_MIN, R_MAX].
    """
    _validate_distance(bond_distance_angstrom)
    return {key: float(spline(bond_distance_angstrom))
            for key, spline in _SPLINES.items()}


def h2_hamiltonian(bond_distance_angstrom: float) -> SparsePauliOp:
    """
    Construct the 2-qubit BK-reduced H2 Hamiltonian as a SparsePauliOp.

    The Hamiltonian form (O'Malley et al., PRX 6, 031007, 2016):

        H(R) = g_I(R)    * II
             + g_Z0(R)   * ZI
             + g_Z1(R)   * IZ
             + g_Z0Z1(R) * ZZ
             + g_Y0Y1(R) * YY
             + g_X0X1(R) * XX

    Qiskit SparsePauliOp string convention (little-endian):
        'ZI' = Z on qubit-1, I on qubit-0 (rightmost)
        'IZ' = I on qubit-1, Z on qubit-0

    Parameters
    ----------
    bond_distance_angstrom : float
        H–H internuclear distance in Ångströms.
        Must lie within [0.20, 2.85] Å.

    Returns
    -------
    SparsePauliOp
        2-qubit Hamiltonian operator (num_qubits == 2).

    Raises
    ------
    ValueError
        If bond_distance_angstrom is outside [R_MIN, R_MAX].
    """
    c = interpolated_coefficients(bond_distance_angstrom)
    pauli_list = [
        (pauli, c[coefficient])
        for pauli, coefficient in zip(PAULI_ORDER, COEFFICIENT_NAMES)
    ]
    return SparsePauliOp.from_list(pauli_list)


def exact_ground_state_energy(bond_distance_angstrom: float) -> float:
    """
    Compute the exact ground-state energy of the H2 Hamiltonian at a
    given bond distance via full matrix diagonalisation.

    Uses the SAME SparsePauliOp returned by h2_hamiltonian() as the
    input Hamiltonian; the exact-reference and VQE paths therefore share
    an identical operator at each bond distance.

    The Hermitian 4×4 matrix is diagonalised with scipy.linalg.eigh
    (symmetric/Hermitian eigensolver), and the smallest eigenvalue is
    returned.

    Parameters
    ----------
    bond_distance_angstrom : float
        H–H internuclear distance in Ångströms.
        Must lie within [0.20, 2.85] Å.

    Returns
    -------
    float
        Minimum eigenvalue of H(R) in Hartree.

    Raises
    ------
    ValueError
        If bond_distance_angstrom is outside [R_MIN, R_MAX].
    """
    op = h2_hamiltonian(bond_distance_angstrom)
    matrix = op.to_matrix()                # 4×4 Hermitian numpy array
    eigenvalues, _ = eigh(matrix)          # sorted ascending, real
    return float(eigenvalues[0])


# ---------------------------------------------------------------------------
# Module-level constants for external use
# ---------------------------------------------------------------------------
PAULI_TERMS: tuple = PAULI_ORDER
"""Ordered Pauli string labels used in the Hamiltonian."""

SOURCE_CITATION: str = (
    "O'Malley, P. J. J., et al., "
    "\"Scalable Quantum Simulation of Molecular Energies,\" "
    "Phys. Rev. X 6, 031007 (2016). arXiv:1512.06860v2."
)
