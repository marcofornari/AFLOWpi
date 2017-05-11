import AFLOWpi

# Create AFLOWpi session
session = AFLOWpi.prep.init('Phonon', 'Perovskite_3',
                            config='./phonon.config')
# Generate a calculation set from a reference input file
allvars={"_AFLOWPI_A_":("Ba","Pb",),
         "_AFLOWPI_B_":("Ti","Zr",),
         "_AFLOWPI_C_":("O",),
         }
calcs = session.scfs(allvars,'phonon.ref',)
#do a vc-relax
calcs.vcrelax()
# change the thresholds in the input files
# changing input for a step are called after
# the call for the step. the change_calcs
# below affect the vc-relax above
calcs.change_input("&control","etot_conv_thr","1.0D-5")
calcs.change_input("&control","forc_conv_thr","1.0D-4")
calcs.change_input("&electrons","conv_thr","1.0D-16")
#do another relax just to be safe
calcs.vcrelax()
#do phonon calculations with 2x2x2 supercell
calcs.phonon(mult_jobs=True,nrx1=2,nrx2=2,nrx3=2,innx=3,
             de=0.003,LOTO=True,field_strength=0.003,
             disp_sym=True,atom_sym=False,)
#plot the phonons
calcs.plot.phonon(postfix='inCM',THz=False)
calcs.plot.phonon(postfix='THz')
# submit the calcs to run
calcs.submit()



