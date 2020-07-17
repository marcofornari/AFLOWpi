import AFLOWpi

# start the AFLOWpirame session
session = AFLOWpi.prep.init('mHT_acbn0', 'zincblende',
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
# calculate the the DOS and PDOS for each 
calcs.dos()
calcs.plot.opdos()
# calculate the bands for each
calcs.bands()
# do the plot the Electronic Band Structure
# and atom projected DOS for each
calcs.plot.bands(DOSPlot='APDOS')
# run the calculation workflow
calcs.acbn0(relax='scf')
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



