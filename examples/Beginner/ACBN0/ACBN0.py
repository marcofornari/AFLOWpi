import AFLOWpi

# start the AFLOWpirame session
session = AFLOWpi.prep.init('ACBN0', 'Perovskite',
                            config='./ACBN0.config')
# choose the values for the keywords in the ref file
allvars{'_AFLOWPI_A_':('Ge',),}
# form the calculation set from ref input and allvars dict
calcs = session.scfs(allvars,'ACBN0.ref')
# relax the structure
calcs.vcrelax()
# calculate the the DOS and PDOS for SrTiO3
calcs.dos()
calcs.plot.opdos(yLim=[-10,10],postfix='without_acbn0')
# calculate the bands for SrTiO3
calcs.bands()
# do the plot the Electronic Band Structure
# and atom projected DOS for SrTiO3
calcs.plot.bands(yLim=[-10,10],DOSPlot='APDOS',
                 postfix='without_acbn0')
# run the ACBN0 pseudo-hybrid functional to
# self-consistently get Hubbard U
calcs.acbn0(thresh=0.1,relax='vc-relax',
            kp_mult=1.5,mixing=0.7)
# calculate the the DOS and PDOS for PBE+U SrTiO3
calcs.dos()
# do the plot the Oribital Proj. DOS for PBE+U SrTiO3
calcs.plot.opdos(yLim=[-10,10],postfix='with_acbn0')
# calculate the bands for PBE+U SrTiO3
calcs.bands()
# do the plot the Electronic Band Structure
# and atom projected DOS for PBE+U SrTiO3
calcs.plot.bands(yLim=[-10,10],DOSPlot='APDOS',
                 postfix='with_acbn0')
# run the calculation workflow
calcs.submit()



