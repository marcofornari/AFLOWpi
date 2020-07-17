import AFLOWpi

# init the frame
session = AFLOWpi.prep.init('Transport', 'Si',
                            config='./transport.config')
# load an already made input file (no AFLOWpi keywords)
calcs = session.from_file('transport.in')
# relax the structure
calcs.vcrelax()
# run calcs for spin polarized NiO for optical and
# transport properties with WanT at 300K and 400K 

tb = calcs.tight_binding(ne=201)
tb.transport(t_min=300,t_max=700,t_step=100,
             en_range=[-2,2])
# plot optical and transport
# properties at 300K and 400K
tb.optical()
tb.plot.transport()
tb.plot.optical()
# submit the workflow to run
calcs.submit()



