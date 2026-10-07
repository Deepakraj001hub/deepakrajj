"""Single-distance statevector VQE for the O'Malley H2 Hamiltonian."""

from dataclasses import dataclass

import numpy as np
from qiskit.circuit.library import efficient_su2
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA, OptimizerResult
from scipy.optimize import minimize as scipy_minimize

from h2_hamiltonian import exact_ground_state_energy, h2_hamiltonian

_MAX_EVALUATIONS = 2000
_OPTIMIZER_TOLERANCE = 1e-5
_INITIAL_TRUST_REGION_RADIUS = 1.0


class _COBYLAWithTermination(COBYLA):
    """COBYLA adapter preserving SciPy termination fields for VQE results."""

    def minimize(self, fun, x0, jac=None, bounds=None) -> OptimizerResult:
        scipy_result = scipy_minimize(
            fun=fun,
            x0=x0,
            method="COBYLA",
            options={
                "maxiter": _MAX_EVALUATIONS,
                "disp": False,
                "rhobeg": _INITIAL_TRUST_REGION_RADIUS,
            },
            tol=_OPTIMIZER_TOLERANCE,
        )
        result = OptimizerResult()
        result.x = scipy_result.x
        result.fun = scipy_result.fun
        result.nfev = scipy_result.nfev
        result.njev = scipy_result.get("njev")
        result.nit = scipy_result.get("nit")
        result.success = bool(scipy_result.success)
        result.status = int(scipy_result.status)
        result.message = str(scipy_result.message)
        return result


@dataclass(frozen=True)
class H2VQEResult:
    """VQE result and exact reference for one H-H bond distance."""

    bond_distance_angstrom: float
    vqe_energy: float
    exact_energy: float
    absolute_error: float
    optimizer: str
    optimizer_configuration: str
    optimizer_evaluations: int
    optimizer_max_evaluations: int
    optimizer_success: bool
    optimizer_status: int
    optimizer_message: str
    maximum_evaluations_reached: bool


def run_h2_vqe(bond_distance_angstrom: float) -> H2VQEResult:
    """Run deterministic statevector VQE for one H2 bond distance.

    Uses a two-qubit EfficientSU2 ansatz with one repetition, a zero-valued
    initial parameter vector, and COBYLA. The exact reference is obtained
    from the existing ``exact_ground_state_energy`` implementation.
    """
    hamiltonian = h2_hamiltonian(bond_distance_angstrom)
    ansatz = efficient_su2(
        2,
        su2_gates=["ry", "rz"],
        entanglement="linear",
        reps=1,
    )
    optimizer = _COBYLAWithTermination(
        maxiter=_MAX_EVALUATIONS,
        tol=_OPTIMIZER_TOLERANCE,
    )
    vqe = VQE(
        estimator=StatevectorEstimator(seed=0),
        ansatz=ansatz,
        optimizer=optimizer,
        initial_point=np.zeros(ansatz.num_parameters),
    )
    vqe_result = vqe.compute_minimum_eigenvalue(hamiltonian)

    vqe_energy = float(np.real(vqe_result.optimal_value))
    exact_energy = exact_ground_state_energy(bond_distance_angstrom)
    return H2VQEResult(
        bond_distance_angstrom=float(bond_distance_angstrom),
        vqe_energy=vqe_energy,
        exact_energy=exact_energy,
        absolute_error=abs(vqe_energy - exact_energy),
        optimizer="COBYLA",
        optimizer_configuration=(
            "maxiter=2000, tol=1e-5, rhobeg=1.0"
        ),
        optimizer_evaluations=int(vqe_result.cost_function_evals),
        optimizer_max_evaluations=_MAX_EVALUATIONS,
        optimizer_success=bool(vqe_result.optimizer_result.success),
        optimizer_status=int(vqe_result.optimizer_result.status),
        optimizer_message=str(vqe_result.optimizer_result.message),
        maximum_evaluations_reached=(
            int(vqe_result.cost_function_evals) >= _MAX_EVALUATIONS
        ),
    )
