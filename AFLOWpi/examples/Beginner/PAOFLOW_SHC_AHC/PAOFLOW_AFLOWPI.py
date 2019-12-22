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

tb.dos()
tb.bands(band_topology=True,fermi_surface=True)

#if t_tensor is omitted (default), all components are calculated
tb.transport(t_tensor=None)

# components for epsilon
d_tensor = [[0,1], # xy component
            [0,2], # xz component
            [2,2]] # zz component

tb.optical(d_tensor=d_tensor)

#use the same as for epsilon
tb.ahc(a_tensor=d_tensor)

# components for SHC
s_tensor=[[0,1,1], # xy component, y spin component
          [0,1,2], # xy component, z spin component
          [2,2,0]] # zz component, x spin component

tb.shc(s_tensor=s_tensor,spin_texture=True)


# submit the workflow to run
calcs.submit()



