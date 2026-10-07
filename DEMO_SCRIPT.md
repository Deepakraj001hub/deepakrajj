# H2 Molecule Energy Estimator — Demo Script

## Slide 1 — Problem → Solution (about 35–40 seconds)

Challenge I5 asks a practical question: how might a pharmaceutical or
materials company estimate a molecule's ground-state energy? We built a
focused demonstration using H2. Its chemistry is represented by a
two-qubit Hamiltonian, and Qiskit's VQE searches for a low-energy state on
a statevector simulator. COBYLA updates the variational parameters. We
compare the resulting estimate with an independent exact classical
eigensolver. In short: molecule, Hamiltonian, VQE, then an exact benchmark.

## Slide 2 — Results: VQE vs Exact (about 50–55 seconds)

This figure shows the energy across 54 bond distances from the published
table. At each point, we ran the same VQE setup and recorded its result,
optimizer termination, and comparison with the exact reference. COBYLA
terminated successfully at 38 points; 16 reached the configured evaluation
limit. That is a 70.4 percent success rate, so convergence was not uniform
across the range. The minimum exact energy is -1.1455991241236425 Hartree
at 0.75 Angstrom. The minimum VQE energy is -1.1455991230750247 Hartree,
also at 0.75 Angstrom. The largest absolute error is 0.00680139457299378
Hartree at 0.25 Angstrom. Red Xs mark unsuccessful optimizer points rather
than hiding them.

## Slide 3 — Business Value → Limits → Next Step (about 25–30 seconds)

This prototype links a molecular Hamiltonian, VQE estimate, exact classical
benchmark, and visible error and convergence analysis. H2 is small enough
for exact classical solution, and simulator results do not establish
quantum advantage. Larger molecules bring chemistry, mapping, circuit-depth,
noise, and scalability challenges. Next, improve convergence handling, then
evaluate LiH; it has not been implemented.
