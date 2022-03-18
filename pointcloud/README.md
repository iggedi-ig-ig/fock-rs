# Pointcloud renderer

This is a pointcloud renderer for the molecular orbitals. This only exists because realtime volumetric rendering isn't feasible.
It works through generating a certain amount of random points around the Molecule and then evaluating the wave function there.
The points are colored based on the sign of the corresponding wave function value, and displayed if their corresponding probability is higher than the current threshold value.

# Example: (all occupied molecular orbitals of Water)
<img alt="plot" height="400" src="https://i.imgur.com/v99VwOd.png" width="530"/>
<img alt="plot" height="400" src="https://i.imgur.com/t5WXtQh.png" width="530"/>
<img alt="plot" height="400" src="https://i.imgur.com/Hp9zfOl.png" width="530"/>
<img alt="plot" height="400" src="https://i.imgur.com/bBeXtjr.png" width="530"/>
<img alt="plot" height="400" src="https://i.imgur.com/ZXKnpe1.png" width="530"/>