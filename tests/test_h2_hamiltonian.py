"""
Unit Tests for src/h2_hamiltonian.py
======================================
Qiskit Fall Fest 2026 – Challenge I5: Molecule Energy Estimator

Tests verify:
  A. Hamiltonian object type is SparsePauliOp.
  B. Hamiltonian has exactly 2 qubits.
  C. Distances inside [0.20, 2.85] Å work without error.
  D. Distances outside [0.20, 2.85] Å raise ValueError.
  E. The canonical published table has the expected shape and range.
  F. Hamiltonian contains exactly the expected six Pauli terms.
  G. exact_ground_state_energy() returns a finite float.
  H. Exact energy is derived from the same Hamiltonian object.
  I. Two published reference points are checked against explicit anchors.
"""

import sys
import os
import math
import unittest
import numpy as np

# Ensure src/ is importable when running from the project root.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from h2_hamiltonian import (
    h2_hamiltonian,
    exact_ground_state_energy,
    interpolated_coefficients,
    R_MIN,
    R_MAX,
    PAULI_TERMS,
)
from h2_table_i_data import (
    COEFFICIENT_NAMES,
    PAULI_ORDER,
    SOURCE_COLUMN_NAMES,
    TABLE_I,
)
from qiskit.quantum_info import SparsePauliOp


# Tolerance for coefficient comparison (round-trip through cubic spline
# at published table nodes should be exact to within floating-point noise).
_COEFF_ABS_TOL = 1e-4   # 0.1 mHartree – sufficient for 4-decimal source data


class TestCanonicalTable(unittest.TestCase):
    def test_table_shape_and_range(self):
        self.assertEqual(len(TABLE_I), 54)
        distances = [row[0] for row in TABLE_I]
        self.assertTrue(all(left < right for left, right in zip(distances, distances[1:])))
        self.assertEqual(distances[0], 0.20)
        self.assertEqual(distances[-1], 2.85)
        self.assertNotIn(0.734, distances)
        self.assertTrue(all(len(row) == 7 for row in TABLE_I))

    def test_declared_pauli_order(self):
        expected = ("II", "ZI", "IZ", "ZZ", "YY", "XX")
        source_order = ("R", "identity", "Z0", "Z1", "Z0Z1", "X0X1", "Y0Y1")
        self.assertEqual(SOURCE_COLUMN_NAMES, source_order)
        self.assertEqual(PAULI_ORDER, expected)
        self.assertEqual(PAULI_TERMS, expected)
        self.assertEqual(len(COEFFICIENT_NAMES), len(PAULI_ORDER))

    def test_all_table_nodes_are_spline_interpolation_nodes(self):
        for row in TABLE_I:
            distance = row[0]
            expected_coefficients = (
                row[1], row[2], row[3], row[4], row[6], row[5]
            )
            with self.subTest(R=distance):
                actual = interpolated_coefficients(distance)
                self.assertEqual(tuple(actual), COEFFICIENT_NAMES)
                for name, expected in zip(COEFFICIENT_NAMES, expected_coefficients):
                    self.assertAlmostEqual(
                        actual[name], expected, delta=_COEFF_ABS_TOL,
                        msg=f"{name} mismatch at R={distance}",
                    )


class TestH2HamiltonianType(unittest.TestCase):
    """A. Object type is SparsePauliOp."""

    def test_returns_sparse_pauli_op(self):
        op = h2_hamiltonian(0.735)
        self.assertIsInstance(
            op, SparsePauliOp,
            "h2_hamiltonian() must return a SparsePauliOp instance."
        )


class TestH2HamiltonianQubits(unittest.TestCase):
    """B. Hamiltonian has exactly 2 qubits."""

    def test_num_qubits_is_2(self):
        op = h2_hamiltonian(0.735)
        self.assertEqual(
            op.num_qubits, 2,
            "H2 Hamiltonian must be defined on exactly 2 qubits."
        )


class TestH2HamiltonianValidDistances(unittest.TestCase):
    """C. Distances inside [R_MIN, R_MAX] work without error."""

    def test_at_r_min(self):
        """Lower boundary R = 0.20 Å must succeed."""
        op = h2_hamiltonian(R_MIN)
        self.assertIsInstance(op, SparsePauliOp)

    def test_at_r_max(self):
        """Upper boundary R = 2.85 Å must succeed."""
        op = h2_hamiltonian(R_MAX)
        self.assertIsInstance(op, SparsePauliOp)

    def test_at_equilibrium(self):
        """Equilibrium distance R ≈ 0.735 Å must succeed."""
        op = h2_hamiltonian(0.735)
        self.assertIsInstance(op, SparsePauliOp)

    def test_mid_range_distances(self):
        """Several mid-range distances must all succeed."""
        for r in [0.30, 0.50, 0.734, 0.75, 1.00, 1.50, 2.00, 2.50, 2.85]:
            with self.subTest(R=r):
                op = h2_hamiltonian(r)
                self.assertIsInstance(op, SparsePauliOp)


class TestH2HamiltonianInvalidDistances(unittest.TestCase):
    """D. Distances outside [R_MIN, R_MAX] raise ValueError."""

    def test_below_r_min(self):
        with self.assertRaises(ValueError):
            h2_hamiltonian(0.10)

    def test_above_r_max(self):
        with self.assertRaises(ValueError):
            h2_hamiltonian(3.00)

    def test_zero(self):
        with self.assertRaises(ValueError):
            h2_hamiltonian(0.0)

    def test_negative(self):
        with self.assertRaises(ValueError):
            h2_hamiltonian(-1.0)

    def test_error_message_contains_range(self):
        """ValueError message must mention the supported range."""
        try:
            h2_hamiltonian(0.05)
            self.fail("Expected ValueError was not raised.")
        except ValueError as exc:
            msg = str(exc)
            self.assertIn("0.20", msg)
            self.assertIn("2.85", msg)


class TestH2HamiltonianCoefficients(unittest.TestCase):
    """E. The Hamiltonian preserves the canonical coefficient/Pauli order."""

    def test_hamiltonian_coefficient_order_matches_canonical_table(self):
        distance = 0.75
        source_row = next(row for row in TABLE_I if row[0] == distance)
        expected_coefficients = dict(zip(COEFFICIENT_NAMES, (
            source_row[1], source_row[2], source_row[3], source_row[4],
            source_row[6], source_row[5],
        )))
        pauli_terms = h2_hamiltonian(distance).to_list()
        self.assertEqual(tuple(label for label, _ in pauli_terms), PAULI_ORDER)
        for (label, actual), expected_label, coefficient_name in zip(
            pauli_terms, PAULI_ORDER, COEFFICIENT_NAMES
        ):
            with self.subTest(pauli=label):
                self.assertEqual(label, expected_label)
                self.assertAlmostEqual(
                    actual, expected_coefficients[coefficient_name],
                    delta=_COEFF_ABS_TOL,
                )


class TestH2HamiltonianPauliTerms(unittest.TestCase):
    """F. Hamiltonian contains exactly the expected six Pauli terms."""

    def test_pauli_labels(self):
        """SparsePauliOp must contain exactly the six expected Pauli strings."""
        op = h2_hamiltonian(0.735)
        labels = set(op.paulis.to_labels())
        expected = set(PAULI_TERMS)
        self.assertEqual(
            labels, expected,
            f"Pauli labels mismatch.\n  Got:      {sorted(labels)}\n"
            f"  Expected: {sorted(expected)}"
        )

    def test_exactly_six_terms(self):
        """SparsePauliOp must have exactly six terms."""
        op = h2_hamiltonian(0.735)
        self.assertEqual(
            len(op), 6,
            f"Expected 6 Pauli terms, got {len(op)}."
        )


class TestExactGroundStateEnergy(unittest.TestCase):
    """G & H. exact_ground_state_energy() returns a finite float from the
    same Hamiltonian."""

    def test_returns_finite_float(self):
        """G. Result must be a finite float."""
        E = exact_ground_state_energy(0.735)
        self.assertIsInstance(E, float)
        self.assertTrue(math.isfinite(E), "Exact energy must be finite.")

    def test_energy_is_negative_at_equilibrium(self):
        """Ground state energy of H2 near equilibrium should be negative."""
        E = exact_ground_state_energy(0.735)
        self.assertLess(E, 0.0, "Ground state energy near equilibrium should be negative.")

    def test_energy_consistent_with_hamiltonian(self):
        """H. Energy is the minimum eigenvalue of h2_hamiltonian()."""
        R = 0.735
        op = h2_hamiltonian(R)
        matrix = op.to_matrix()
        eigenvalues = np.linalg.eigvalsh(matrix)
        E_ref = float(np.min(eigenvalues))
        E_func = exact_ground_state_energy(R)
        self.assertAlmostEqual(
            E_func, E_ref, places=10,
            msg="exact_ground_state_energy() must equal the minimum eigenvalue "
                "of h2_hamiltonian()."
        )

    def test_energy_at_multiple_distances(self):
        """Energy must be finite and real at several distances."""
        for r in [0.20, 0.735, 1.00, 1.50, 2.00, 2.50]:
            with self.subTest(R=r):
                E = exact_ground_state_energy(r)
                self.assertTrue(math.isfinite(E))

    def test_energy_decreases_from_compressed_to_equilibrium(self):
        """Energy should decrease from compressed (R=0.20) toward equilibrium."""
        E_compressed = exact_ground_state_energy(0.20)
        E_equilibrium = exact_ground_state_energy(0.735)
        self.assertGreater(
            E_compressed, E_equilibrium,
            "Energy at compression (R=0.20) should exceed equilibrium energy."
        )


class TestPublishedReferencePoint(unittest.TestCase):
    """I. Check two explicit published anchors against the canonical table."""

    _ANCHORS = {
        0.20: (2.8489, 0.5678, -1.4508, 0.6799, 0.0791, 0.0791),
        0.75: (0.2252, 0.3435, -0.4347, 0.5716, 0.0910, 0.0910),
    }

    def test_published_anchor_values(self):
        canonical_by_distance = {row[0]: row[1:] for row in TABLE_I}
        for distance, published_values in self._ANCHORS.items():
            with self.subTest(R=distance):
                self.assertEqual(canonical_by_distance[distance], published_values)
                actual = interpolated_coefficients(distance)
                internal_values = (
                    published_values[0], published_values[1],
                    published_values[2], published_values[3],
                    published_values[5], published_values[4],
                )
                for name, expected in zip(COEFFICIENT_NAMES, internal_values):
                    self.assertAlmostEqual(
                        actual[name], expected, delta=_COEFF_ABS_TOL,
                        msg=f"{name} mismatch at published R={distance}",
                    )


if __name__ == "__main__":
    unittest.main(verbosity=2)
