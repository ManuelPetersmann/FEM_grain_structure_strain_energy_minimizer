# FEM grain structure strain energy minimizer

Given a grain structure, where each grain can have 1...N grain different grain strains.
<img src="https://raw.githubusercontent.com/ManuelPetersmann/FEM_grain_structure_strain_energy_minimizer/master/img17.png" width="400">
<!-- ![](https://raw.githubusercontent.com/ManuelPetersmann/FEM_grain_structure_strain_energy_minimizer/master/img17.png | width=100)  -->

and under suitable boundary conditions (e.g. periodic boundary conditions, or embedding the grain structure into a matrix like shown here)
<img src="https://raw.githubusercontent.com/ManuelPetersmann/FEM_grain_structure_strain_energy_minimizer/master/img16.png" width="400">
<!-- ![test](https://raw.githubusercontent.com/ManuelPetersmann/FEM_grain_structure_strain_energy_minimizer/master/img16.png | width=100)  -->

the given phython scripts transform each grain, one by one, trying out all possible transformation strains a grain can have in the following energy-minimizing manner:

<img src="https://raw.githubusercontent.com/ManuelPetersmann/FEM_grain_structure_strain_energy_minimizer/master/procedure_IEMA.JPG" width="400">

## Reference
If you've found this useful please consider referencing this repository in your own work
```
@article{Petersmann_2017,
	doi = {10.1088/1361-651x/aa5ab4},
	url = {https://doi.org/10.1088/1361-651x/aa5ab4},
	year = 2017,
	month = {feb},
	publisher = {{IOP} Publishing},
	volume = {25},
	number = {3},
	pages = {035004},
	author = {M Petersmann and T Antretter and T Waitz and F D Fischer},
	title = {A new approach predicting the evolution of laminated nanostructures{\textemdash}martensite in {NiTi} as an example},
	journal = {Modelling and Simulation in Materials Science and Engineering},
	abstract = {A model for laminated nanostructures, combining classical energy minimization with full-field finite element calculations in a computationally fully automated manner, is set up and used to quantitatively analyse the interaction of grains via self-accommodation of their transformation strains. The well established B2–B19’ martensitic phase transformation in nanocrystalline NiTi is treated as an exemplary case to demonstrate our new framework. A systematic search for an optimal energy minimizing transformation path is employed within a full-field model, including crystallographic transformation strains and fully anisotropic elastic constants, by using the Python scripting language. The microstructure is updated based on previous calculation results. The underlying incremental free energy minimization criterion naturally reproduces the transformation kinetics. The sequence of grains subjected to transformation as well as the selection of martensitic variants within the grains are obtained yielding the evolution of the total interface energy as well as the strain energy, dominating our approach.}
}
```
