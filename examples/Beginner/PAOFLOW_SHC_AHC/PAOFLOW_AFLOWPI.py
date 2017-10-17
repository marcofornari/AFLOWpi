import AFLOWpi

# init the frame
session = AFLOWpi.prep.init('AHC_SHC', 'fe',
                            config='./PAOFLOW_AFLOWPI.config')
# load an already made input file (no AFLOWpi keywords)
calcs = session.from_file('PAOFLOW_PW.in')
# relax the structure
calcs.scf()
# run calcs for spin polarized NiO for optical and
# transport properties with WanT at 300K and 400K 

tb = calcs.tight_binding()
tb.bands(topology=True)
# plot optical and transport
# properties at 300K and 400K
tb.ahc()
tb.shc()


# submit the workflow to run
calcs.submit()



