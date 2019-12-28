import AFLOWpi

# start the AFLOWpirame session
session = AFLOWpi.prep.init('mHT', 'zincblende',
                            config='./mHT.config')
# choose the values for the keywords in the ref file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Al',),
_AFLOWPI_B_ = ('As',),)
# form the calculation set from ref input and allvars dict
calcs = session.from_file(['mHT.in'])
# relax the structure
calcs.vcrelax()
#updated cell parameters
# calculate the the DOS and PDOS for each 
calcs.vcrelax()
calcs.vcrelax()
calcs.dos()
calcs.plot.opdos()
# calculate the bands for each
calcs.bands()
# do the plot the Electronic Band Structure
# and atom projected DOS for each
calcs.plot.bands(DOSPlot='APDOS')
# run the calculation workflow
#calcs.submit()



