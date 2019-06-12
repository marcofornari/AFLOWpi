
import AFLOWpi

# Define the values to the keywords
# in the reference input file

# Create AFLOWpi session
session = AFLOWpi.prep.init('testing', 'GaAs',config='./epsilon.config')

# Generate a calculation set from a reference input file
calcs = session.from_file("GaAs.in")



calcs.change_input("K_POINTS","__content__","8 8 8 1 1 1")

calcs.scf()
calcs.epsilon(offdiag=True,occ=True,jdos=True)


calcs.submit()



