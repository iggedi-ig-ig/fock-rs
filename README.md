# fock-rs
Restricted Hartree Fock in rust

This is not fully working, as for example d-orbitals don't work correctly (the orbitals / energies that are formed don't make any sense)
Also the scf-routine returns slightly lower values for energies than expected (which doesn't make any sense)

Features:
- pointcloud renderer for molecular orbitals
- custom boys function implementation which is decently fast with good enough accuracy
- parsing of the BSE json format of basis sets using serde-rs
