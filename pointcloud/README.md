This is a pointcloud renderer for the molecular orbitals. This only exists because realtime volumetric rendering isn't feasible.
It works through generating a certain amount of random points around the Molecule and then evaluating the wave function there.
The points are colored based on the sign of the corresponding wave function value, and displayed if their corresponding probability is higher than the current threshold value.

Example (pi orbital of benzene):
![plot](https://i.imgur.com/OLop4OA.png)

Note that this example is from an earlier state, so don't mind the incorrect energy values.