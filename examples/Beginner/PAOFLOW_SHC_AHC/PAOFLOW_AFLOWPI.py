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
# plot optical and transport
# properties at 300K and 400K
tb.transport()
tb.optical()
tb.dos()
tb.bands(band_topology=True,fermi_surface=True)
tb.ahc()

s_tensor=[[0,1,1],[0,1,2],[2,2,0]]
tb.shc(s_tensor,spin_texture=True)


# submit the workflow to run
calcs.submit()



