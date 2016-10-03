import AFLOWpi
# start the AFLOWpirame session
session = AFLOWpi.prep.init('mHT_TB', 'zincblende',config='./mHT_TB.config')
# choose the values for the keywords in the ref file
allvars={}
allvars.update(
_AFLOWPI_A_ = ('Ga','In'),
_AFLOWPI_B_ = ('As','P'),)
# form the calculation set from ref input and allvars dict
calcs = session.scfs(allvars,'mHT_TB.ref')
# relax the structure
calcs.vcrelax()
#calculate the PAO-TB Hamiltonian
tb_ham = calcs.tight_binding()
#calculate the DOS with PAO-TB
tb_ham.dos()
tb_ham.plot.opdos()
#calculate the electronic band structure with PAO-TB
tb_ham.bands()
tb_ham.plot.bands(DOSPlot='APDOS')
#submit the calcs to run
calcs.submit()



