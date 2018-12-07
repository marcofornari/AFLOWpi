import AFLOWpi

# Define the values to the keywords
# in the reference input file

# Create AFLOWpi session
session = AFLOWpi.prep.init('testing', 'cDFT_Fe3O4',config='./cDFT_test.config')


# Generate a calculation set from a reference input file
calcs = session.from_file(["Fe3O4.in"])

ch1=6.80644
ch2=5.95592

oxid_states={'Fe1':ch1,'Fe2':ch2}

calcs.force_oxidation(oxid_states,conv_thr=0.01)

calcs.submit()

