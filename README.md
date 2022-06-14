# Treball de Fi de Grau de Física
##### Universitat Autònoma de Barcelona
## Studying Weak Turbulence in non-linear Wave Dynamics

#### Miquel de la Iglesia Martí
##### Supervisat per Dr. Carlos Fernández Sopuerta, Institut de Ciències de l’Espai (ICE, CSIC i IEEC)
There are two folders: One for the one-dimensional simulations and the other for two-dimensional ones. Inside each folder there are three Python files: the corresponding to the main simulation, another one to generate the initial conditions and a third one to plot the Power Spectra that we obtain from the evolution.

For each simulation one has to introduce the equation that you would like to evolve, together with the time duration, the number of spatial grid steps. One can also choose whether or not you would like to generate a .GIF file out of the evolution data.

_IMPORTANT!_<br />
The number of the spatial grid points for the simulation must be the same as the one used to generate the initial conditions. If you would like to change the discretization,  you should first run the code to generate the initial conditions with the new number of spatial grid points. 
