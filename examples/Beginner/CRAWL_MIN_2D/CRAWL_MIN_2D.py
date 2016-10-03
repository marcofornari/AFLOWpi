import AFLOWpi

# Create AFLOWpi session
session = AFLOWpi.prep.init('CRAWL_MIN', 'graphene',
                            config='./CRAWL_MIN_2D.config')
# Generate a calculation set from a list of input files
calcs = session.from_file(['CRAWL_MIN_2D.in',])
# do a crawling minimization of the lattice 
# parameters for HEX lattice. C is fixed.
calcs.crawl_min(constraint=['fixed','c'],
                initial_variance=0.05,
                thresh=0.0001,grid_density=10,
                mult_jobs=False)
#do a vc-relax at the end just for verification
calcs.vcrelax()
# submit the calcs to run
calcs.submit()



