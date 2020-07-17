import AFLOWpi
# start the AFLOWpirame session
session = AFLOWpi.prep.init('AFLOWpi_tests', 
                            'thermal_ZB',config='./thermal.config')
# choose the values for the keywords in the ref file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Al',),
_AFLOWPI_B_ = ('P',),)
# form the calculation set from ref input and allvars dict
calcs = session.scfs(allvars,'thermal.ref',)
# relax the structure
calcs.vcrelax()
# change the thresholds in the input files
# changing input for a step are called after
# the call for the step. the change_calcs
# below affect the vc-relax above
calcs.change_input("&control","etot_conv_thr","1.0D-5")
calcs.change_input("&control","forc_conv_thr","1.0D-4")
calcs.change_input("&electrons","conv_thr","1.0D-12")
calcs.vcrelax()

calcs.phonon(mult_jobs=True,nrx1=2,nrx2=2,nrx3=2,innx=2,
             LOTO=False,disp_sym=True,atom_sym=True,)

# do thermal calc with forward difference finite 
# difference calc of gruneisen parameter
calcs.thermal(delta_volume=0.04,mult_jobs=True,nrx1=2,
              nrx2=2,nrx3=2,innx=2,de=0.01,LOTO=False,
              disp_sym=True,atom_sym=True,central_diff=True)

# FYI: this won't have LOTO splitting
calcs.plot.phonon(postfix='inCM',THz=False,runlocal=False)
#plot k_L vs. Temp
calcs.plot.thermal_cond(temp_range=[250,800])
#plot frequency resolved (projected part still needs work)
calcs.plot.gruneisen()

calcs.submit()



