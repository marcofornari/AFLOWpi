import AFLOWpi

# start the AFLOWpirame session
session = AFLOWpi.prep.init('AFLOW_TEST', 'partial_inputs',
                            config='./partial_in.config')

# form the calculation set from ref input and allvars dict
#calcs = session.from_file(['TET.in','BCT_1.in','BCT_2.in'],
#calcs = session.from_file(['CUB.in','FCC.in','BCC.in'],
calcs = session.from_file(['BCC.in'],
#                          reffile='./common_values_template.ref'
                          reffile='/mnt/home/supka1arCMICH/AFLOWpi/examples/Intermediate/Partial_Input_Files/common_values_template.ref')
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



