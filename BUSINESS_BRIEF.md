# H2 Molecule Energy Estimator
## Industry Challenge I5 — Business Brief

### 1. Customer Problem

A pharmaceutical or materials firm needs a molecule's ground-state energy
to help assess molecular structures and materials. Challenge I5 asks how a
quantum-computing workflow could estimate that energy.

### 2. Proposed Solution

This prototype studies H2 with a two-qubit Hamiltonian and Qiskit VQE,
running on a statevector simulator. COBYLA adjusts the ansatz parameters;
an exact classical eigensolver supplies the reference. The experiment
evaluates both approaches across the published bond-distance range.

### 3. What Was Demonstrated

The saved sweep and its summary report:

- 54 bond-distance points.
- 38 successful optimizer terminations and 16 maximum-evaluation
  terminations: a 70.4% optimizer success rate.
- Minimum exact energy: -1.1455991241236425 Ha at 0.75 Angstrom.
- Minimum VQE energy: -1.1455991230750247 Ha at 0.75 Angstrom.
- Largest absolute error: 0.00680139457299378 Ha at 0.25 Angstrom.

### 4. Classical Benchmark

The independent exact eigensolver is the reference for this small
Hamiltonian. The VQE estimates are compared with that reference to measure
energy error and inspect optimizer termination. No claim is made that the
quantum method is faster or better than the classical solution.

### 5. Business Value

The prototype demonstrates a workflow:

**molecular Hamiltonian -> VQE estimate -> classical benchmark ->
error/convergence analysis**

This workflow can help study the behavior of a quantum estimation pipeline
before considering larger chemistry problems. No ROI, cost savings,
runtime improvement, or production capability is established.

### 6. Limitations and Risks

- H2 is a very small demonstration system, and an exact classical solution
  is feasible for it.
- 16 of 54 VQE runs reached the configured evaluation limit; convergence
  is not uniform across the bond-distance range.
- Simulator results do not establish quantum advantage.
- Larger molecules introduce additional chemistry, mapping, circuit-depth,
  noise, and scalability challenges not addressed by this demonstration.

### 7. Recommended Next Step

Improve optimizer and convergence handling, then evaluate a somewhat larger
chemistry example such as LiH while retaining exact-classical benchmarking.
LiH has not been implemented in this project.

### 8. Executive Takeaway

- **Built:** an H2, two-qubit Qiskit VQE workflow using statevector
  simulation and an exact classical reference.
- **Demonstrated:** a 54-point bond-distance experiment with recorded
  energies, errors, and optimizer termination outcomes.
- **Unproven:** reliable convergence across the full range, scalability to
  larger chemistry problems, and any quantum advantage.
