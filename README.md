**Installation instructions for AFLOWpi**

**PREREQUISITES:**
Python     >= 3.5 
NumPy      >= 1.8.0
SciPy      >= 0.14.0
matplotlib >= 1.3.0

The following free scientific Python distributions include these prerequisites:

*Continuum Anaconda*
https://www.anaconda.com/products/individual

*Enthought Canopy*
https://assets.enthought.com/downloads/edm/

**For ACBN0 and PAOFLOW capabilities:**
PAOFLOW >= 2.0.10
mpi4py >= 2.0

**INSTALLATION:**
To install AFLOWpi, execute this command in the AFLOW$\pi$ source directory:

1. Clone this repository and enter the cloned directory
2. python setup.py install --user

with the --user flag optional depending on the setup of your python distribution.

**AFLOW$\pi$ can also be installed easily onto linux based systems via pip:**
pip install AFLOWpi

**GETTING STARTED:**
 The examples directory in the AFLOWpi source directory contains example scripts. 
- AFLOW$\pi$ config must be modified to point to the proper executable directories of the calculation engine(s) being used. 
- Once the AFLOWpi config file is properly setup for the computer environment, to start the example, run them as python executables. 


**FOR ACBN0:**
Quantum Espresso must be modified and recompiled for ACBN0 to work properly. The files needed are located in engine_mods directory. 
