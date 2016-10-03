import AFLOWpi

# Define the values to the keywords
# in the reference input file

# Create AFLOWpi session
session = AFLOWpi.prep.init('CRAWL_MIN', 'perovskite',
                            config='./CRAWL_MIN_TET.config')

# Generate a calculation set from a list of input files
calcs = session.from_file(['CRAWL_MIN_TET.in',])
#do a crawling minimization of the lattice parameters for TET
calcs.crawl_min(initial_variance=0.1,thresh=0.0001,
                grid_density=10,mult_jobs=True)
#do a vc-relax at the end just for verification
calcs.vcrelax()
# submit the calcs to run
calcs.submit()



