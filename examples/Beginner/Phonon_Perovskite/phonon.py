import AFLOWpi
import glob

# Create AFLOWpi session
session = AFLOWpi.prep.init('Phonon', 'Perovskite',
                            config='./phonon.config')

# get a list of input files using glob
in_files = glob.glob("input_files/*.in")

# make a calculation set from those files
calcs = session.from_file(in_files,
                          filename_as_dirname=True)
#do a vc-relax
calcs.vcrelax()
# change the thresholds in the input files
# changing input for a step are called after
# the call for the step. the change_calcs
# below affect the vc-relax above
calcs.change_input("&control","etot_conv_thr","1.0D-5")
calcs.change_input("&control","forc_conv_thr","1.0D-4")
calcs.change_input("&electrons","conv_thr","1.0D-12")
#do another relax just to be safe
calcs.vcrelax()
#do phonon calculations with 2x2x2 supercell
calcs.phonon(mult_jobs=True,nrx1=2,nrx2=2,nrx3=2,innx=2,
             de=0.003,LOTO=True,field_strength=0.003,
             field_cycles=4,disp_sym=True,atom_sym=True,)
#plot the phonons
calcs.plot.phonon(postfix='inCM',THz=False)
calcs.plot.phonon(postfix='THz')
# submit the calcs to run
calcs.submit()



