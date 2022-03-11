This contains the main self consistent field logic using [nalgebra](https://nalgebra.org/). 
This is probably the crate that contains the error that is causing energies to be off.
My best guesses are:
- the multi electron tensor is using symmetry wrongly
- I'm using nalgebras SymmetricEigen wrongly and therefore get eigenvector / value pairs that are slightly off
- Some matrix that should be symmetric / hermitian is not for some reason
