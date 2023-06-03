# fock-rs
Restricted Hartree Fock in rust

This is just for educational purposes as only little thought went into optimizing everything. Especially 
the electron tensor is really inefficient, as I don't employ any fancy optimizations like integral screening, etc.

# Features:
- pointcloud and volumetric renderer for molecular orbitals
- custom boys function implementation, which is fast enough but pretty inaccurate
- parsing of the BSE json format of basis sets using serde-rs

# Problems:
- the scf cycle is not perfectly stable. It converges for smaller molecules with simple basis sets, but struggles for larger systems
- the boys function implementation is pretty inaccurate for higher order orbitals, which is why basis sets with polarization can be unstable
- poor performance / very high impact of ERIs
