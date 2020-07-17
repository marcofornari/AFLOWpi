import AFLOWpi
# start the AFLOWpirame session
session = AFLOWpi.prep.init('Workflow_Append',SET='Si',
                            config='./do_calcs.config')
# choose the values for the keywords in the ref file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Si',),
_AFLOWPI_B_ = ('Si',),)
# form the calculation set from ref input and allvars dict
calcs = session.scfs(allvars,'do_calcs.ref')

# relax the structure
calcs.vcrelax()
calcs.crawl_min(initial_variance=0.10,thresh=0.001,
                grid_density=10,mult_jobs=True)
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
# do thermal properties. if phonons previously calculated 
# and structure/charge density don't change then thermal
# properties will skip initial phonon calculation at V_0
#NOTE: this will give bad results. it's just for testing
calcs.thermal(mult_jobs=True,nrx1=1,nrx2=1,nrx3=1,innx=1,
              LOTO=False,disp_sym=True)
#do elastic constants
calcs.elastic()
#submit the calculations to run.
calcs.submit()
