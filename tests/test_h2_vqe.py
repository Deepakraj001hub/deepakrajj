"""Focused tests for the single-distance H2 statevector VQE."""

import math
import os
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import h2_vqe
from h2_hamiltonian import exact_ground_state_energy


class TestH2VQE(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with patch(
            "h2_vqe.exact_ground_state_energy",
            wraps=exact_ground_state_energy,
        ) as exact_solver:
            cls.results = {
                distance: h2_vqe.run_h2_vqe(distance)
                for distance in (0.75, 0.50, 1.50)
            }
        cls.exact_solver_calls = exact_solver.call_args_list
        cls.result = cls.results[0.75]

    def test_module_imports_and_result_structure(self):
        self.assertTrue(callable(h2_vqe.run_h2_vqe))
        self.assertEqual(
            set(self.result.__dataclass_fields__),
            {
                "bond_distance_angstrom",
                "vqe_energy",
                "exact_energy",
                "absolute_error",
                "optimizer",
                "optimizer_configuration",
                "optimizer_evaluations",
                "optimizer_max_evaluations",
                "optimizer_success",
                "optimizer_status",
                "optimizer_message",
                "maximum_evaluations_reached",
            },
        )

    def test_vqe_and_reference_energies_are_finite(self):
        for result in self.results.values():
            with self.subTest(R=result.bond_distance_angstrom):
                self.assertTrue(math.isfinite(result.vqe_energy))
                self.assertTrue(math.isfinite(result.exact_energy))

    def test_exact_energy_comes_from_existing_solver(self):
        self.assertEqual(len(self.exact_solver_calls), 3)
        self.assertEqual(
            [call.args for call in self.exact_solver_calls],
            [(0.75,), (0.50,), (1.50,)],
        )
        for distance in (0.75, 0.50, 1.50):
            with self.subTest(R=distance):
                self.assertEqual(
                    self.results[distance].exact_energy,
                    exact_ground_state_energy(distance),
                )

    def test_absolute_error_is_computed_from_returned_energies(self):
        for result in self.results.values():
            with self.subTest(R=result.bond_distance_angstrom):
                self.assertEqual(
                    result.absolute_error,
                    abs(result.vqe_energy - result.exact_energy),
                )

    def test_three_representative_distances_include_termination_metadata(self):
        self.assertEqual(set(self.results), {0.50, 0.75, 1.50})
        for distance, result in self.results.items():
            with self.subTest(R=distance):
                self.assertEqual(result.bond_distance_angstrom, distance)
                self.assertEqual(result.optimizer, "COBYLA")
                self.assertEqual(
                    result.optimizer_configuration,
                    "maxiter=2000, tol=1e-5, rhobeg=1.0",
                )
                self.assertGreater(result.optimizer_evaluations, 0)
                self.assertEqual(result.optimizer_max_evaluations, 2000)
                self.assertIsInstance(result.optimizer_success, bool)
                self.assertIsInstance(result.optimizer_status, int)
                self.assertTrue(result.optimizer_message)
                self.assertTrue(result.optimizer_success)
                self.assertFalse(result.maximum_evaluations_reached)
                self.assertEqual(
                    result.maximum_evaluations_reached,
                    result.optimizer_evaluations >= result.optimizer_max_evaluations,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
