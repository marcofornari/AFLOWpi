
import AFLOWpi

# Define the values to the keywords
# in the reference input file

# Create AFLOWpi session
session = AFLOWpi.prep.init('testing', 'GaAs',config='./everything.config')

# Generate a calculation set from a reference input file
calcs = session.from_file("GaAs.in")

# relax the structure and prepare for Finite Diff. Phonons

calcs.change_input("K_POINTS","__content__","4 4 4 1 1 1")
calcs.vcrelax()
calcs.vcrelax()
calcs.change_input("&electrons","conv_thr","1.e-10")
calcs.vcrelax()

# calculate bands/pdos without acbn0
calcs.dos()
calcs.bands(nk=200)
calcs.plot.opdos()
calcs.plot.bands(DOSPlot='APDOS',)


# run acbn0 to get hubbard U
calcs.acbn0(relax='vc-relax',thresh=0.01)

calcs.vcrelax()
calcs.vcrelax()

# calculate bands/pdos/transport with ACBN0
calcs.dos()
calcs.bands(nk=200)
calcs.plot.opdos(yLim=[-4,4])
calcs.plot.bands(DOSPlot='APDOS',yLim=[-4,4])


calcs_tb=calcs.tight_binding(proj_thr=0.9,
                             smearing='m-p',
                             kp_factor=2.0,
                             tb_kp_factor=8.0,
                             emin=-1.0,
                             emax=1.0)

calcs_tb.dos()
calcs_tb.bands()
calcs_tb.transport(t_min=300,
                   t_max=700,
                   t_step=100)
calcs_tb.optical()

calcs_tb.plot.optical(x_range=[-4,4],postfix="testing")
calcs_tb.plot.transport(x_range=[-4,4],postfix="testing")
calcs_tb.plot.opdos(yLim=[-4,4],postfix="testing")
calcs_tb.plot.bands(DOSPlot="APDOS",yLim=[-4,4],postfix="testing")

calcs.change_input("K_POINTS","__content__","4 4 4 1 1 1")
# do phonon and elastic constants
calcs.phonon(LOTO=False,mult_jobs=True,raman=False,disp_sym=True,de=0.005)
calcs.plot.phonon(postfix='in_THz',THz=True)
calcs.plot.phonon(postfix='before_thermal_call',THz=False)


#calcs.elastic(mult_jobs=False)

calcs.thermal(mult_jobs=True,disp_sym=True,LOTO=False,central_diff=True,de=0.005)

calcs.plot.gruneisen()
calcs.plot.phonon(postfix='after_thermal_call',THz=True)
calcs.plot.thermal_cond(temp_range=[300,700])


calcs.change_input("K_POINTS","__content__","4 4 4 1 1 1")
# prep inputs for LSDA calcs
calcs.change_input("&system","starting_magnetization(1)","-0.5")
calcs.change_input("&system","starting_magnetization(2)","0.5")
calcs.change_input("&system","nspin",'2')    
calcs.change_input("&system",'smearing','"mv"')  
calcs.change_input("&system",'occupations','"smearing"')  
calcs.change_input("&system",'degauss','0.001')  

calcs_tb=calcs.tight_binding(proj_thr=0.9,
                             smearing='m-p',
                             kp_factor=2.0,
                             tb_kp_factor=8.0,
                             emin=-1.0,
                             emax=1.0)

calcs_tb.dos()
calcs_tb.bands(nk=2000)
calcs_tb.transport(t_min=300,
                   t_max=700,
                   t_step=100)



calcs_tb.plot.transport(x_range=[-4,4],postfix="noncolin")
calcs_tb.plot.opdos(yLim=[-4,4],postfix="noncolin")
calcs_tb.plot.bands(DOSPlot="APDOS",yLim=[-4,4],postfix="noncolin")

calcs.dos()
calcs.bands(nk=200)
calcs.plot.opdos(yLim=[-4,4])
calcs.plot.bands(DOSPlot='APDOS',yLim=[-4,4])

#prep inputs for noncolin calcs
calcs.change_pseudos('./pseudo_rel')
calcs.change_input("&system","lda_plus_u_kind",'1')
calcs.change_input("&system","lspinorb",".true.")
calcs.change_input("&system","noncolin",".true.")
calcs.change_input("&system","nspin",None) 
calcs.change_input("&system",'smearing',None)  
calcs.change_input("&system",'occupations',None)  
calcs.change_input("&system",'degauss',None)  

calcs_tb=calcs.tight_binding(proj_thr=0.9,
                             smearing='m-p',
                             kp_factor=2.0,
                             tb_kp_factor=8.0,
                             emin=-1.0,
                             emax=1.0)

calcs_tb.dos()
calcs_tb.bands(nk=2000)
calcs_tb.transport(t_min=300,
                   t_max=700,
                   t_step=100)

calcs_tb.ahc()
calcs_tb.shc()


calcs_tb.plot.transport(x_range=[-4,4],postfix="noncolin")
calcs_tb.plot.opdos(yLim=[-4,4],postfix="noncolin")
calcs_tb.plot.bands(DOSPlot="APDOS",yLim=[-4,4],postfix="noncolin")
calcs_tb.plot.shc()
calcs_tb.plot.ahc()

calcs.dos()
calcs.bands(nk=1000)
calcs.plot.opdos(yLim=[-4,4])
calcs.plot.bands(DOSPlot='APDOS',yLim=[-4,4])

calcs.shake_atoms(dist=0.1,weight=True)
calcs.relax()

calcs.submit()



