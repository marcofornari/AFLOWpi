import AFLOWpi

# init the frame for the three sets
# NOTE: Config must be the same for all
session_cub = AFLOWpi.prep.init('CRT','CUB'
                        ,config='./CRTs.config')
session_tet = AFLOWpi.prep.init('CRT','TET'
                        ,config='./CRTs.config')
session_rho = AFLOWpi.prep.init('CRT','RHO'
                        ,config='./CRTs.config')
# it's time to tell what you want to do:
# all tuples of strings or tuple of numbers 
allvars={'_AFLOWPI_KK_': ('6 6 6',),
	 '_AFLOWPI_KS_': ('1 1 1',),
	 '_AFLOWPI_A_' : ('Sr','Ba',),
	 '_AFLOWPI_B_' : ('Zr','Ti',),
	 '_AFLOWPI_C_' : ('O',),}
# prepare the dictionary of dictionaries
# with all the calculations
calcs_cub = session_cub.scfs(allvars, './5CUB.ref')
calcs_tet = session_tet.scfs(allvars, './5TET.ref')
calcs_rho = session_rho.scfs(allvars, './5RHO.ref')
#relax the calcs for step 1 before ACBNO
calcs_cub.vcrelax()
calcs_tet.vcrelax()
calcs_rho.vcrelax()
# submit the calcs to run
calcs_cub.submit()
calcs_tet.submit()
calcs_rho.submit()
