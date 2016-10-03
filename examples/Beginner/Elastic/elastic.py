import AFLOWpi

# Define the values to the keywords
# in the reference input file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Si',),
_AFLOWPI_B_ = ('Si',),)
# Create AFLOWpi session
session = AFLOWpi.prep.init('Elastic', 'Si',
                            config='./elastic.config')
# Generate a calculation set from a reference input file
calcs = session.scfs(allvars,'elastic.ref')
#calcs.scfuj(relax='vc-relax')
calcs.vcrelax()
# do the elastic constants with the ElaStic Package 
# Install: http://exciting-code.org/elastic
calcs.elastic(mult_jobs=False,num_dist=10,eta_max=0.001)
# submit the calcs to run
calcs.submit()



