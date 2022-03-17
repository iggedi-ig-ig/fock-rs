# fock-rs
Restricted Hartree Fock in rust

This is just for educational purposes as only little thought went into optimizing everything. Especially 
the electron tensor is really inefficient, as I don't employ any fancy optimizations like integral screening, etc.

# Features:
- pointcloud renderer for molecular orbitals
- custom boys function implementation which is decently fast with good enough accuracy
- parsing of the BSE json format of basis sets using serde-rs

# Problems:
- the scf cycle currently gets stuck in a cycle pretty often, at least for more complex systems like benzene. This could probably be countered by using a better initial guess for the density matrix
- poor performance / very high impact of ERIs