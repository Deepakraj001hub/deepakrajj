# H2 Molecule Ground-State Energy Estimator — Qiskit VQE

## Industry problem

The Qiskit Fall Fest Industry Track I5 challenge asks how a pharmaceutical
or materials company can estimate a molecule's ground-state energy, a
quantity relevant to comparing molecular structures and screening candidate
materials. This project demonstrates that workflow for the small hydrogen
molecule (H2); it is a prototype, not a production chemistry platform.

## Approach

The experiment uses a two-qubit Hamiltonian for H2 and Qiskit's Variational
Quantum Eigensolver (VQE) with a statevector simulator and the COBYLA
classical optimizer. An independent exact classical eigensolver provides a
reference energy. The bond-distance experiment evaluates both methods at
each published distance in the source table.

The VQE uses a two-qubit `efficient_su2` ansatz with `ry`/`rz` gates, linear
entanglement, one repetition, and a zero-valued initial parameter vector.
COBYLA uses `maxiter=2000`, `tol=1e-5`, and `rhobeg=1.0`.

## Data provenance

The canonical Hamiltonian data in `src/h2_table_i_data.py` preserves the
published values from Table I of O'Malley et al., *Scalable Quantum
Simulation of Molecular Energies*, Physical Review X 6, 031007 (2016),
arXiv:1512.06860v2, DOI: 10.1103/PhysRevX.6.031007. The table contains 54
published bond distances from 0.20 Angstrom through 2.85 Angstrom at
0.05 Angstrom spacing.
Published table rows are kept as tabulated; intermediate Hamiltonian
queries use spline interpolation and are not additional published data.

The paper describes a minimal basis of Hartree-Fock orbitals. This README
does not assign an unsupported basis-set label.

## Method

For each bond distance, the workflow is:

1. Build the two-qubit Hamiltonian for H2.
2. Prepare the parameterized ansatz.
3. Evaluate its energy with the statevector estimator.
4. Use classical COBYLA optimization to vary the ansatz parameters.
5. Record the optimized VQE energy.
6. Compare it with the independent exact classical eigensolver and compute
   the absolute error.

## Results

The completed sweep contains 54 points:

| Measure | Result |
| --- | ---: |
| Successful optimizer terminations | 38 / 54 (70.4%) |
| Unsuccessful / maximum-evaluation terminations | 16 / 54 |
| Minimum exact energy | -1.1455991241236425 Ha at 0.75 Angstrom |
| Minimum VQE energy | -1.1455991230750247 Ha at 0.75 Angstrom |
| Largest absolute error | 0.00680139457299378 Ha at 0.25 Angstrom |

The numerical summary is in `results/h2_vqe_summary.json`, and the
point-by-point results and optimizer metadata are in
`results/h2_vqe_sweep.csv`.

## Figures

- [Energy vs. bond distance](figures/h2_energy_vs_bond_distance.png)
- [VQE absolute error](figures/h2_vqe_absolute_error.png)

Unsuccessful optimizer points are explicitly marked in both figures and
remain included in the plotted data.

## Limitations

- 16 of 54 points reached the configured evaluation limit; convergence is
  not uniform across the bond-distance range.
- H2 is a very small demonstration system.
- These simulator results do not demonstrate quantum advantage.
- An exact classical eigensolver is feasible for this small Hamiltonian.
- This prototype demonstrates a workflow; it is not evidence of
  production-scale quantum advantage.

## Reproducibility

The verified environment uses Python 3.13.14 and the pinned dependencies in
`requirements.txt` (including Qiskit 2.5.2 and qiskit-algorithms 0.4.0).
From the repository root, create and activate the environment and install
the pinned requirements:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

The saved sweep results and figures are included in the repository; setup
and testing do not require rerunning the VQE sweep.

## Repository structure

- `src/` — canonical published Hamiltonian data, Hamiltonian construction,
  VQE implementation, and sweep runner.
- `tests/` — unit tests for the Hamiltonian and VQE result contract.
- `results/` — validated sweep CSV and JSON summary.
- `figures/` — energy and absolute-error plots from the sweep CSV.

## Hackathon framing

**Predict -> Build -> Benchmark -> Explain**

This project predicts H2 ground-state energies with VQE, builds the
Hamiltonian and ansatz workflow, benchmarks against an exact classical
reference, and explains optimizer outcomes and limitations alongside the
results.
