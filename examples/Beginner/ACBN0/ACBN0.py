import AFLOWpi

# start the AFLOWpirame session
session = AFLOWpi.prep.init('ACBN0', 'Si2',
                            config='./ACBN0.config')


# form the calculation set from ref input and allvars dict
calcs = session.from_file('Si2.in')
# relax the structure
calcs.vcrelax()
# calculate the the DOS and PDOS for Si2
calcs.dos()
calcs.plot.opdos(en_range=[-10,10],postfix='without_acbn0')
# calculate the bands for Si2
calcs.bands(nk=200)
# do the plot the Electronic Band Structure
# and atom projected DOS for Si2
calcs.plot.bands(en_range=[-10,10],DOSPlot='APDOS',
                 postfix='without_acbn0')
# run the ACBN0 pseudo-hybrid functional to
# self-consistently get Hubbard U
calcs.acbn0(thresh=0.1,relax='scf',kp_factor=2.0)            
# calculate the the DOS and PDOS for PBE+U Si2
calcs.vcrelax()
calcs.dos()
# do the plot the Oribital Proj. DOS for PBE+U Si2
calcs.plot.opdos(en_range=[-10,10],postfix='with_acbn0')
# calculate the bands for PBE+U Si2
calcs.bands(nk=200)
# do the plot the Electronic Band Structure
# and atom projected DOS for PBE+U Si2
calcs.plot.bands(en_range=[-10,10],DOSPlot='APDOS',
                 postfix='with_acbn0')
# run the calculation workflow
calcs.submit()



