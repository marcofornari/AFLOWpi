import AFLOWpi

#init AFLOWpirame
session = AFLOWpi.prep.init('pseudotesting','zincblende',
                            config='./pseudotesting.config')
# Choose 3 species for the test
allvars={'_AFLOWPI_A_':['Ge','Si','C',]}
#pseudotesting ref contains zincblende FCC structure
calcs=session.scfs(allvars,"pseudotesting.ref")
# run the pseudotests on the 3 species in the
# zincblende structure initial variance of A is
# +-10% with 4 varations tested per iteration of
# the crawling minimization (more variations gives
# better fit)
calcs.pseudo_test_brute([15,16,17,18,20,25,35,50,60,70],
                        sampling=['10 10 10 1 1 1'],
                        initial_variance=0.10,
                        conv_thresh=0.01,
                        min_thresh=0.01,grid_density=4,
                        mult_jobs=False)
# vc-relax at the converged cutoff and sampling 
# to see how much the structure will change vs. 
# brute force optimization
calcs.vcrelax()
#submit the calcs to run
calcs.submit()
