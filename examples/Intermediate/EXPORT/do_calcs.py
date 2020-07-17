import AFLOWpi
# start the AFLOWpirame session
session = AFLOWpi.prep.init('export_example', 'zincblende_mHT',config='./do_calcs.config')
# choose the values for the keywords in the ref file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Ga',),
_AFLOWPI_B_ = ('As',),)
# form the calculation set from ref input and allvars dict
calcs = session.scfs(allvars,'do_calcs.ref')

# relax the structure
calcs.vcrelax()
#do acbn0 to get the hubbard U vals
calcs.acbn0(relax='vc-relax')
# relax the structure
calcs.vcrelax()
calcs.crawl_min(initial_variance=0.10,thresh=0.001,grid_density=10,mult_jobs=True)
#calculate the PAO-TB Hamiltonian
tb_ham = calcs.tight_binding()
#calculate the DOS with PAO-TB
tb_ham.dos()
#calculate the electronic band structure with PAO-TB
tb_ham.bands()
tb_ham.plot.opdos()
tb_ham.plot.bands()
tb_ham.effmass()
tb_ham.plot.bands(DOSPlot='APDOS')
#do thermal properties 
calcs.thermal()
#do elastic constants
calcs.elastic()
#do dos with bloch representation
calcs.dos()
calcs.plot.dos()
calcs.plot.opdos()
#do bands with bloch representation
calcs.bands()
calcs.plot.bands()
calcs.plot.bands(DOSPlot='DOS')
calcs.plot.bands(DOSPlot='APDOS')
#do phonon
calcs.phonon(mult_jobs=True,LOTO=True,field_strength=0.001)
calcs.plot.phonon()
#do transport properties
tb_ham.transport(temperature=[300,400])
tb_ham.plot.transport()
tb_ham.optical()
tb_ham.plot.optical()
#submit the calcs to run
calcs.submit()
