import AFLOWpi

# Define the values to the keywords
# in the reference input file

# Create AFLOWpi session
session = AFLOWpi.prep.init('testing', 'environ',config='./environ_test.config')


# Generate a calculation set from a reference input file
calcs = session.from_file(["si.in"])


#not really needed..just to test
calcs.scf()

# doesnt do much right now..
calcs.environ()



calcs.submit()

