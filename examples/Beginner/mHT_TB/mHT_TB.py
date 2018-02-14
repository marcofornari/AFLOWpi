import AFLOWpi

# init the frame
session = AFLOWpi.prep.init('mHT_TB', 'ZB',
                            config='./mHT_TB.config')

allvars={}
allvars.update(
_AFLOWPI_A_ = ('Al','Ga'),
_AFLOWPI_B_ = ('As','P'),)

calcs = session.scfs(allvars,'mHT_TB.ref')

# relax the structure
calcs.vcrelax()
calcs.update_cell()

calcs.vcrelax()
calcs.update_cell()

# run calcs for spin polarized NiO for optical and
# transport properties with WanT at 300K

# adaptive smearing is incompatable with calc of kappa and seebeck
# proj_thr set low for HT..0.9 is liwest probably for Colusite HT
# nscf_kp_grid = kp_factor*(scf_kp_grid)
# PAOTB_kp_grid = tb_kp_factor*(scf_kp_grid)
# I suggest tb_kp_factor=12 for Colusite HT (16+ if you want smoother transport plots)
tb = calcs.tight_binding(smearing=None,proj_thr=0.85,kp_factor=2.0,tb_kp_factor=8.0)

# dos+pdos
tb.dos(projected=True)

# do tb bands
tb.bands(band_topology=False,fermi_surface=False)

#if t_tensor is omitted (default), all components are calculated
tb.transport(t_tensor=None)

# components for epsilon
# must have xx,yy,zz for plotting epsilon
d_tensor = [[0,1], # xy component
            [0,2], # xz component
            [0,0], # xx component
            [1,1], # yy component
            [2,2]] # zz component

tb.optical(d_tensor=d_tensor)

tb.plot.bands(yLim=[-1.5,1.5],DOSPlot="APDOS")
tb.plot.opdos(yLim=[-1.5,1.5])

# plot sigma,kappa,seebeck
# x_range is mu range
tb.plot.transport(x_range=[-1.5,1.5])

# plot epsilon
tb.plot.optical(x_range=[-1.5,1.5])

# submit the workflow to run
calcs.submit()



