
import AFLOWpi
# load logs from each set of calcs
calcs_CUB=AFLOWpi.prep.init('CRT','CUB',config='./CRTs.config').load(1)
calcs_TET=AFLOWpi.prep.init('CRT','TET',config='./CRTs.config').load(1)
calcs_RHO=AFLOWpi.prep.init('CRT','RHO',config='./CRTs.config').load(1)


# plot the deltaE between the one of the sets of calculations and the rest of
# them. In this case we are plotting the energy difference between the cubic
# (nondistorted) lattice and the (distored) rhombohedral and tetragonal
# lattices. 'titleArray' are the titles for the color bars for each of the 
# two distortions. (first distortion being rhombohedral and second tetragonal)
AFLOWpi.plot.grid_plot([calcs_CUB,calcs_TET,calcs_RHO],'_AFLOWPI_A_','_AFLOWPI_B_',
                   zaxis_title=['Tet $\Delta$E Ry','Rho $\Delta$E Ry',],
                   plot_title='$ABO_{3}$ Distortion $\Delta$E',yAxisStr='B',
                   xAxisStr='A',)

