# VULCAN
#### Authors: Shang-Min (Shami) Tsai ####
Photochemical kinetics for exoplanetary atmospheres, a fast and easy-to-use python code.
VULCAN is implemented with the equilibrium chemistry code [FastChem](https://github.com/exoclime/FastChem) created by Daniel Kitzmann, Joachim Stock, to initialise a state in chemical equilibrium. 

The theory paper of the first version of VULCAN (without photochemistry) can be found here: https://arxiv.org/abs/1607.00409

This is currectly a release candidate version. Any questions or feedbacks is welcome and can be sent to [Shami Tsai](mailto:shang-min.tsai@physics.ox.ac.uk)

## Requirements
VULCAN is developed with Python 3 but has been tested compatible with Python 2.7. It is advised to run it on Python 3 if possible.
Two very useful tools to set up python environments:
[Pip](https://pip.pypa.io/en/stable/) - package installer for Python
[Anaconda](https://docs.continuum.io/) - virtual environment manager

VULCAN requires the following python packages:
- numpy
- scipy
- Sympy
- matplotlib
- PIL/Pillow (optional: for interactive plotting)
and the embeded [FastChem](https://github.com/exoclime/FastChem) requires a standard C++ compiler, such as g++ or Clang.

If any of the python packages are missing, you can install the full SciPy Stack via Pip, e.g.
```bash
python -m pip install --upgrade pip
``` 
```bash
pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
```
The above commands update pip and install SciPy via pip. Further information can be found at http://www.scipy.org/install.html

PIL or Pillow is a plotting library. If installed, the plots will be conveniently shown by the os-built-in image viewer. See https://github.com/python-pillow/Pillow for more information.  

## Quick Demo

Let's dive in and see how chemical kinetics works!

First, go to the ```/fastchem_vulcan``` folder to compile [FastChem](https://github.com/exoclime/FastChem)(equilibrium chemistry code) by running 
```
make
```

After compiling finished, go back to the main directory of VULCAN and run 
```
python vulcan.py -np 
```

You should see the default model for HD 189733b starts running with real-time plotting. This will take about 10-15 minutes to complete depending on your comuputer. Why not make some coffee, sit back and enjoy watching the chemistry in action!

After the run finished with steady state reached, we can plot the results from the output (stored in ```/output``` by default). Run the  plotting script ```plot_vulcan.py``` in the folder ```plot_vpy```, followed by three arguments: **{output path} {comma-separated species} {plotname}**. For example,
```
python plot_vulcan.py output/HD189.vul H2O,CH4,CO,CO2,NH3,HCN hd189
```
will plot the output file "HD189.vul" for the chosen species: H2O,CH4,CO,CO2,NH3,HCN and save the plot named "hd189" in the ```/output``` folder. 


Now you may want to try a different T-P input or elemental abundances. All these settings are prescreibed in ```vulcan_cfg.py```. For example, find and edit
```python
const_Kzz = 1.E7
```
and 
```python
C_H = 6.0618E-4
``` 
for a smaller diffusion (K<sub>zz</sub>) and carbon rich (C/O=1) run. Set use_live_plot = False if you wish to switch off the real-time plotting. More detailed instruction can be found in the following sections. Have fun!

## Full instruction

### Structure
```
├── VULCAN/
│   ├── atm/
│   ├── /thermo/
│   │   ├──/NASA9/
│   │   ├── all_compose.txt
│   │   ├── gibbs_text.txt
│   ├── build_atm.py  
│   ├── chem_funs
│   ├── op.py
│   ├── prepipe.py
│   ├── plot_vulcan.py
│   ├── phy_const.py
│   ├── store.py
│   ├── vulcan.py
│   ├── vulcan_cfg.py
``` 

`/atm/` : storing input atmospheric files  
`/thermo/` : storing chemical kinetics networks and thermodynamic data
   
`/NASA9/` : storing the NASA-9 polynomials for the Gibbs free energy of every species

`/thermo/all_compose.txt` : basic physical property e.g. number of atoms and molecular weight
`build_atm.py` : modules to construct the atmospheric structure and to read in reaction rates  
`chem_funs.py` : the functions of chemical sources, Jacobian matrix and the equilibrium constants    
`CHO_netowrk.txt` : the default C-H-O kinetics network  
`op.py` : modules for the main computation  
`prepipe.py` : pre-pipeline routine to produce `chem_funs.py`    
`plot_vulcan.py` : read-in and plotting script  
`phy_const.py` : physical constants  
`store.py` : modules to store all the variables  
`vulcan.py` : the main file of VULCAN  
`vulcan_cfg.py` : the configuration file for VULCAN  

```vulcan_cfg.py``` includes all the settings and parameters, e.g. the atmospheric parameters, the elemental abundance etc. Typically this is the only file you need to edit for each specific run. A summary of every setting is listed in ```vulcan_cfg_readme.txt```. 

### Input Files
The key input files of VULCAN include the chemical network and the atmospheric file.
```CHO_network.txt``` is the deafult reaction network including carbon, hydrogen, and oxygen species. It is constructed for  thermochemistry from 500 to 2500 K using a reduced network with 29 gaseous species and less than 300 reactions to benefit for efficiency. The rate coefficients A, B, C are written in A, B, C as in the Arrhenius formula k = A T^B exp(-C/T) (See the section of editing or using a different chemical network.)  The input temperature-pressure profiles are placed in the `/atm` folder by default. As the form of the included example files for HD 189733b and HD 209458b, the first line is always commented for units, and the second line specifies the column names: **Pressure	Temp** or **Pressure	Temp  	Kzz** (if Kzz-P profile is provided). So the file consists of two columns without K<sub>zz</sub> and three columns with K<sub>zz</sub>. When K<sub>zz</sub> is provided, set
```python
Kzz_prof = 'file'
```
in the configuration file ```vulcan_cfg.py```.
  
### Making chem_funs.py based on the assigned chemical network
VULCAN is developed in a flexible way that the chemical network is _not_ hard coded. Instead, ```make_chem_funs.py``` generates all the required funtions from the input chemical network (e.g. ```NCHO_photo_netowrk.txt```) into ```chem_funs.py```. By default, ```make_chem_funs.py``` is always run prior to the main code. This step (which takes a few seconds) can be skiped by adding the agument ```-n```while running vulcan in the command line:
```
python vulcan.py -n 
```
only when change networks

It is important to not skipping this step after making a change of the chemical network. You can examine all the readable information, like the list of reactions and species in ```chem_funs.py```, being updated.

### Editing or Using a different chemical network
You can edit the default netowrk, to remove or add reactions, to change rate constats, etc. You can also use a different chemical network, as long as it is in the same format as ```CHO_netowrk.txt```. Noted that however, changing or using a different chemical network is not warranted, unrealistic values could lead to numerical issues. In the network file, the reactions should be writen in the form of [ A + B -> C + D ], including the square brackets. Only the forward reactions are listed, since VULCAN reverses the forward reactions to obtain the reverse reactions using the thermodynamic data. The reaction number, i.e. **id**, is irrelevent as it will be automatically generated (and writing into the network file) while running ```prepipe.py```. Three-body or dissociation reactions should be separately listed after the comment line 
```python

### Boundary Condistions ###


# 3-body and Disscoiation Reactions
```
Next, make sure all the species are included in the ```NASA9``` folder. If not, they need to be add manually by looking over ```nasa9_2002_E.txt``` and save them into txt files with the same name as used in the network. The format of the NASA 9 polynomials is as follows
```
a1 a2 a3 a4 a5
a6 a7 0. a8 a9
```
Noted that a7 and a8 are separated by 0. The first two rows are for low temperature (200 - 1000 K) and the last two rows are for high temperature (1000 - 6000 K).

Finally, run ```prepipe.py``` as described before.

## License
VULCAN is distributed under the terms of the GNU General Public License (GPL) license. For more information, see ```GPL_license.txt``` in the main directory.

## Remarks
The project is financially support from the Center for Space and Habitability (CSH), the PlanetS NCCR framework and the Swiss- based MERAC Foundation.
The Exoclime Simulation Platform ([ESP][1]) develops a set of open-source codes
for research on exoplanets. The three parts of the ESP are
  - [HELIOS][2] radiative transfer and retrieval,
  - [THOR][3] atmospheric fluid dynamics,
  - [VULCAN][4] atmospheric chemistry.

[1]: http://www.exoclime.net
[2]: https://github.com/exoclime/HELIOS
[3]: https://github.com/exoclime/THOR
[4]: https://github.com/exoclime/VULCAN
