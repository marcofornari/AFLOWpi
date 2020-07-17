import AFLOWpi

# initiate AFLOWpi
session = AFLOWpi.prep.init('Workflow_Append', 
                            SET='Si',
                            config='./do_calcs.config')

# load from the last step in the previously completed workflow
calcs=session.load(8)

#start appending more steps to the end of the workflow
calcs.dos()
calcs.plot.dos()
calcs.plot.opdos()
#do bands with bloch representation
calcs.bands()
calcs.plot.bands()
calcs.plot.bands(DOSPlot='DOS')
calcs.plot.bands(DOSPlot='APDOS')
#do transport properties
calcs.transport(temperature=[300,400])
calcs.plot.transport()
#submit the calcs to run
calcs.submit()
