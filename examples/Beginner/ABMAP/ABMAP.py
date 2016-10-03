import AFLOWpi
import numpy

# init a session of AFLOWpirame 
session=AFLOWpi.prep.init('ABMAP',config='./ABMAP.config')
#define the input
BaTiO3_input='''
 &control
    calculation = 'scf',
    prefix = '_AFLOWPI_PREFIX_',
    restart_mode='from_scratch', 
    PSEUDO_DIR ='./',outdir ='./',
 /
 &system
    ibrav=1,celldm(1) = 7.51915, 
    nat=5,ntyp=3, ecutwfc = 40.,
 /
 &electrons
 /
ATOMIC_SPECIES
_AFLOWPI_A_  _AFLOWPI_AMASS_ _AFLOWPI_APSEUDO_ 
_AFLOWPI_B_  _AFLOWPI_BMASS_ _AFLOWPI_BPSEUDO_ 
_AFLOWPI_C_  _AFLOWPI_CMASS_ _AFLOWPI_CPSEUDO_ 
ATOMIC_POSITIONS {crystal}
_AFLOWPI_A_       _AFLOWPI_AX_   0.000000  0.000000
_AFLOWPI_B_       _AFLOWPI_BX_   0.500000  0.500000
_AFLOWPI_C_       0.500000   0.500000  0.000000
_AFLOWPI_C_       0.500000   0.000000  0.500000
_AFLOWPI_C_       0.000000   0.500000  0.500000
K_POINTS automatic
6 6 6 1 1 1'''

allvars={}
# do numpy linspace to get an evenly spaces set of 5 values 
# from 0.5 to 0.7 for B site and 0.0 to 0.2 for A site
allvars.update(
_AFLOWPI_A_ = ('Ba',),
_AFLOWPI_B_ = ('Ti',),
_AFLOWPI_C_ = ('O',),
_AFLOWPI_AX_ =  numpy.linspace(0.0, 0.10, 5),
_AFLOWPI_BX_ =  numpy.linspace(0.5, 0.60, 5),)
# make the 5x5 grid of calculations with 
# A and B atoms shifting in in X axis
calcs=session.scfs(allvars,BaTiO3_input)
# do an scf to get the energy grid
calcs.scf()
# start calculations 
#calcs.submit()

