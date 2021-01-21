**Installation instructions for AFLOWpi**

**PREREQUISITES:**

- Python     >= 3.5 
- NumPy      >= 1.8.0
- SciPy      >= 0.14.0
- matplotlib >= 1.3.0

The following free scientific Python distributions include these prerequisites:

- *Continuum Anaconda*
https://www.anaconda.com/products/individual
- *Enthought Canopy*
https://assets.enthought.com/downloads/edm/

**For ACBN0 and PAOFLOW capabilities:**

- PAOFLOW >= 2.0.10
https://github.com/marcobn/PAOFLOW
- mpi4py >= 2.0

**INSTALLATION:**

To install AFLOWpi, execute this command in the AFLOWpi source directory:

1. Clone this repository and enter the cloned directory
2. python setup.py install --user

with the --user flag optional depending on the setup of your python distribution.

**AFLOWpi can also be installed easily onto linux based systems via pip:**

pip install AFLOWpi

**GETTING STARTED:**
 
The examples directory in the AFLOWpi source directory contains example scripts. 
- AFLOWpi config must be modified to point to the proper executable directories of the calculation engine(s) being used. 
- Once the AFLOWpi config file is properly setup for the computer environment, to start the example, run them as python executables. 


**FOR ACBN0:**

Quantum Espresso must be modified and recompiled for ACBN0 to work properly. The files needed are located in engine_mods directory. 

Use of AFLOWpi should reference:

A.R. Supka, T.E. Lyons, L. Liyanage, P. D'Amico, R. Al Rahal Al Orabi, S. Mahatara, P. Gopal, C. Toher, D. Ceresoli, A. Calzolari, S. Curtarolo, M. Buongiorno Nardelli, and M. Fornari, AFLOW𝛑: A minimalist approach to high-throughput ab initio calculations including the generation of tight-binding hamiltonians, Computational Materials Science, 136 (2017) 76-84. doi:10.1016/j.commatsci.2017.03.055

Use of PAOFLOW should reference:

M. Buongiorno Nardelli, F. T. Cerasoli, M. Costa, S Curtarolo,R. De Gennaro, M. Fornari, L. Liyanage, A. Supka and H. Wang, PAOFLOW: A utility to construct and operate on ab initio Hamiltonians from the Projections of electronic wavefunctions on Atomic Orbital bases, including characterization of topological materials, Comp. Mat. Sci. vol. 143, 462 (2018). doi:10.1016/j.commatsci.2017.11.034



