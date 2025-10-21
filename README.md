# VULCAN
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)\
Photochemical kinetics for exoplanetary atmospheres, a fast and easy-to-use python code.

## Credits
* Shang-Min (Shami) Tsai
* Harrison Nicholls

This distribution of VULCAN contains a number of performance and usability improvements.

The theory papers of VULCAN can be found here: [Tsai et al. 2021](https://arxiv.org/abs/2108.01790) (with photochemistry) and [Tsai et al. 2017](https://arxiv.org/abs/1607.00409) (without photochemistry).
This is currently a release candidate version. Any questions or feedbacks is welcome and can be sent to [Shami Tsai](mailto:shang-min.tsai@physics.ox.ac.uk)

* Demo with realtime plotting:\
![Running with realtime plotting](demo/demo.gif)

## Quick Demo

Let's dive in and see chemical kinetics in action!

First, go to the `fastchem_vulcan/` folder to compile [FastChem](https://github.com/shami-EEG/FastChem) by running
```
make
```

After compilation has finished, go back to the main directory of VULCAN and run
```
python vulcan.py
```

You should see the default model starts running with real-time plotting.
This will take about 10-15 minutes to complete depending on your comuputer.

Now you may want to try a different T-P input, changing the elemental abundances or
vertical mixing. All these settings are set in ```config.py```. For example, find and edit
```python
const_Kzz = 1.E7
```
and
```python
C_H = 6.0618E-4
```
for a weaker vertical mixing (K<sub>zz</sub>) and carbon rich (C/O=1) run. Set use_live_plot = False if you wish to switch off the real-time plotting (why whould you though?). More detailed instruction can be found in the following sections. Have fun!

The object in this `config.py` file can be edited at runtime and passed around as a variable.

## Repository structure

`/atm/`: storing input atmospheric files
`/fastchem_vulcan/`: Fastchem (equilibirum chemistry code) which is used to initialse the compositions
`/output/`: storing the output files
`/plot/`: storing the output plots
`/thermo/`: storing chemical kinetics networks and thermodynamic data
`/thermo/NASA9/`: storing the NASA-9 polynomials for the Gibbs free energy of every species
`/thermo/all_compose.txt`: basic compositional properties e.g. number of atoms and molecular weight
`/thermo/gibbs_text.txt`: a text file used by make_chem_funs.py to generate chem_funs.py
`build_atm.py`: modules to construct the atmospheric structure based on the input and to set up the initial compositions
`chem_funs.py`: the functions of chemical sources/sinks, Jacobian matrix and the equilibrium constants
`op.py`: all the modules for the numerical operations, e.g. computing reaction rates, ODE solvers etc.
`make_chem_funs.py`: the routine that runs first to produce the required `chem_funs.py` based on the assigned chemical network
`phy_const.py`: physical constants
`store.py`: modules to store all the variables
`vulcan.py`: the top-level main script of VULCAN
`config.py`: the configuration file for VULCAN

Typically ```config.py``` is the only file you need to edit for each specific run. If you want to look inside or modify the code, `store.py` is where almost all classes and variables are declared.

## Configuration File
<strong>All the settings and parameters, e.g. the atmospheric parameters, the elemental abundance etc, are prescribed in ```vulcan_cfg.py```</strong>. Typically this is the only file you need to edit for each specific run. A useful cheatsheet describing what every parameter does can be found in ```vulcan_cfg_readme.txt```. The configuration files used for the model validation in [Tsai et al. 2021](https://arxiv.org/abs/2108.01790) are also provided in the cfg_examples folder.

## Input Files
The key input files of VULCAN include the chemical network, atmospheric T-P profile, and stellar flux. ```NCHO_photo_network.txt``` is the deafult reaction network including nitrogen, carbon, hydrogen, and oxygen species. It is validated from ~ 500 to 3000 K with about 60 gaseous species and 700 reactions.
The rate coefficients A, B, C are written in A, B, C as in the Arrhenius formula k = A T^B exp(-C/T).
The input temperature-pressure(-Kzz) profile is required when Kzz_prof is set to 'file' in vulcan_cfg.py and is placed in the `/atm` folder by default. The first line in the T-P file is commented for units, and the second line must specifies the column names: **Pressure	Temp** or **Pressure	Temp  	Kzz** (Kzz is optional). So the file consists of two columns without K<sub>zz</sub> and three columns with K<sub>zz</sub>.
See the included T-P files of HD 189733b and HD 209458b in `/atm` for example.
The stellar UV flux is stored in /atm/stellar_flux, with the first column being weavelength in nm and the second column	being flux in ergs/cm**2/s/nm.
The thermodynamics data and cross sections are stored in /thermo/NASA9 and /thermo/photo_cross, respectively. Change at your own risk!
If constant fluxes for certain species are used, the files are also placed in /atm, in the format of species, flux (cm-2 s-1), and deposite velocity (cm s-1).

## Editing or Using a different chemical network
VULCAN is developed in a flexible way that the chemical network is _not_ hard coded. Instead, ```make_chem_funs.py``` generates all the required funtions from the input chemical network (e.g. ```NCHO_photo_network.txt```) into ```chem_funs.py```.
You can edit the default netowrk, to remove or add reactions, to change rate constats, etc. You can also use a different chemical network, as long as it is in the same format as the defalut ones. That is, the reactions should be writen in the form of [ A + B -> C + D ], including the square brackets.
By default, ```make_chem_funs.py``` is always called prior to the main code to produce ```chem_funs.py``` based on the new chemical network . This step (which takes a few seconds) can be skiped by adding the agument ```-n```while running vulcan in the command line:
```
python vulcan.py -n
```
However, it is important NOT to skipping this step after making a change of the chemical network.

Noted that changing or using a different chemical network is not foolproof -- unrealistic values could lead to numerical issues. You can also see that only the forward reactions are listed, since VULCAN reverses the forward reactions to obtain the reverse reactions using the thermodynamic data.
So next, make sure all the species are included in the ```NASA9``` folder. If not, they need to be added manually by looking over ```nasa9_2002_E.txt``` or ```new_nasa9.txt```, which can also be found in ```NASA9```. Save the coefficients in a text file with the same name as used in the network (e.g. CO2.txt). The format of the NASA 9 polynomials is as follows
```
a1 a2 a3 a4 a5
a6 a7 0. a8 a9
```
Here, a7 and a8 are separated by 0. The first two rows are for low temperature (200 - 1000 K) and the last two rows are for high temperature (1000 - 6000 K).\

The reaction number, i.e. **id**, is irrelevent as it will be automatically generated (and writing into the network file) while calling ```make_chem_funs```. Three-body or dissociation reactions should also be separately listed after the comment line as the default network.
After changing the network, you can examine all the readable information, like the list of reactions and species in ```chem_funs.py```, being updated while running python vulcan.py (without -n argument).

## Boundary Conditions
If both use_topflux and use_botflux in vulcan_cfg.py are set to False, it will use the default boundary condition -- zero flux boundary i.e. nothing coming in or out. When use_topflux = True, it reads the file prescribed in top_BC_flux_file as the incoming/outgoing flux at the top boundary. Similarly, when use_botflux = True, the file prescribed in bot_BC_flux_file is read in for the surface pressure and sinks at the bottom boundary. In addition, you can also use the dictionary use_fix_sp_bot to set fixed mole fraction at the surface. e.g. use_fix_sp_bot = {'CO2': 0.01} sets the surface CO<sub>2</sub> mixing ratio to 0.01.

## Coupling to a self-consistent climate model
VULCAN can be self-consistently coupled to a atmosphere climate model (AGNI) which solves for the atmospheric temperature profile given the composition calculated by VULCAN. AGNI uses correlated-k radiative transfer and mixing-length convection in order to determine realistic TP- and Kzz-profiles. More information on AGNI can be found [here](https://www.h-nicholls.space/AGNI/).

To enable chemistry-climate calculations, install AGNI on your machine and make it available inside the VULCAN folder (e.g. `VULCAN/AGNI/agni.jl`). Then set `agni_call_frq` to a value >0 in the VULCAN config object and run VULCAN like normal.

## Reading Output Files
Run ```plot_vulcan.py``` within ```plot_py```
```
python plot_vulcan.py [vulcan output] [species] [plot name] [-h (for plotting height)]
```
will read vulcan output (.vul files) can plot the species profiles. Species should be sepreated by commas without space. Plot is in pressure by diffcult and can be changed to height by adding "-h".
