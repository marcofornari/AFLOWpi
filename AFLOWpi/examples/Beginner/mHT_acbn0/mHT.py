import AFLOWpi

# start the AFLOWpirame session
session = AFLOWpi.prep.init('mHT', 'zincblende',
                            config='./mHT.config')
# choose the values for the keywords in the ref file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Ga','In'),
_AFLOWPI_B_ = ('As','P'),)
# form the calculation set from ref input and allvars dict
calcs = session.scfs(allvars,'mHT.ref')
# relax the structure
calcs.vcrelax()
calcs.pseudo_test_brute([15,16,17,18,20,25,35,50,60,70],
                        sampling=['10 10 10 1 1 1'],
                        initial_variance=0.10,
                        conv_thresh=0.01,
                        min_thresh=0.01,
                        grid_density=4)
# calculate the the DOS and PDOS for each 
calcs.dos()
calcs.plot.opdos()
# calculate the bands for each
calcs.bands()
# do the plot the Electronic Band Structure
# and atom projected DOS for each
calcs.plot.bands(DOSPlot='APDOS')
# run the calculation workflow
calcs.scfuj(relax='vc-relax')
# relax the structure
calcs.vcrelax()
# calculate the the DOS and PDOS for each 
calcs.dos()
calcs.plot.opdos(postfix='w_acbn0')
# calculate the bands for each
calcs.bands()
# do the plot the Electronic Band Structure
# and atom projected DOS for each
calcs.plot.bands(DOSPlot='APDOS',postfix='w_acbn0')
# submit the calculations
calcs.submit()



