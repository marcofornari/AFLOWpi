import AFLOWpi
import glob

# start the AFLOWpirame session
session = AFLOWpi.prep.init('AFLOW_TEST', 'partial_inputs',
                            config='./partial_in.config')


files=glob.glob("inputs/*.in")

calcs = session.from_file(files,filename_as_dirname=True,
                          reffile='./common_values_template.ref')
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
#calcs.submit()



