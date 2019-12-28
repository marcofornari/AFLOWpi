import AFLOWpi

# Create AFLOWpi session
session = AFLOWpi.prep.init('Berry', 'perovskite',
                            config='./berry.config')
# Generate a calculation set from a reference input file
allvars={"_AFLOWPI_A_":("Ba","Pb",),
         "_AFLOWPI_B_":("Ti","Zr",),
         "_AFLOWPI_C_":("O",),
         }
calcs = session.scfs(allvars,'berry.ref',)



#do a vc-relax
calcs.vcrelax()

#shake atom positions randomly to break centrosymmetry
calcs.shake_atoms(dist=0.1)

#do calc for Berry phase
calcs.berry()

# submit the calcs to run
#calcs.submit()



