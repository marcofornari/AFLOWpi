import AFLOWpi

# Load the logs from the calcs with U
CUB_SESSION =  AFLOWpi.prep.init('ABMAP',
                                 config='./ABMAP.config')
cub_after=CUB_SESSION.load(1)
# grab the energy from output
cub_after=AFLOWpi.retr.grabEnergyOut(cub_after)
# Plot map of Energy as Ti moves in XY plane
AFLOWpi.plot.interpolatePlot(cub_after,'_AFLOWPI_AX_','_AFLOWPI_BX_',
                             zaxis='Energy',fileName='BaTiO3_enMap.pdf',
                             delta_min=False,
                             xaxisTitle="Ti Position: A Axis (crys. coord.)",
                             yaxisTitle="Ti Position: B Axis (crys. coord.)",
                             zaxisTitle="Energy (Ry)",
                             title='Energy of Cubic\nPerovskite BaTiO$_3$, $A=4.029\AA$')
