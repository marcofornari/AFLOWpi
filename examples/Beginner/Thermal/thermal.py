import AFLOWpi

# Define the values to the keywords
# in the reference input file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Si',),
_AFLOWPI_B_ = ('Si',),)
# Create AFLOWpi session
session = AFLOWpi.prep.init('Thermal', 'Si',
                            config='./thermal.config')
# Generate a calculation set from a reference input file
calcs = session.scfs(allvars,'thermal.ref')
# relax the structure and prepare for Finite Diff. Phonons
#calcs.vcrelax()
calcs.vcrelax()
# calculate one phonon frequency
calcs.thermal(mult_jobs=True,nrx1=3,nrx2=3,nrx3=3,innx=2,
              field_strength=0.001,de=0.003,LOTO=True,
              delta_volume=0.03,disp_sym=True)
# plot phonon dispersion and DOS
calcs.plot.phonon()


# submit the calcs to run
calcs.submit()



