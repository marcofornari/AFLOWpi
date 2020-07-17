import AFLOWpi
import glob

# Create AFLOWpi session
session = AFLOWpi.prep.init('AFLOWPI_EXAMPLES', 'HALL',
                            config='./Hall.config')

# get a list of input filenames
files=glob.glob("inputs/*.in")

# Generate a calculation set from pwscf input files
calcs = session.from_file(files,filename_as_dirname=True)

#optimize structure with scalar-relativistic pseudos
calcs.vcrelax()

# generate PAO-TB hamiltonian with w/o SOC
calcs_tb=calcs.tight_binding(proj_thr=0.95,tb_kp_factor=4,
                             symmetrize=False)

#calculate bands and DOS with no SOC
calcs_tb.dos()
calcs_tb.bands()

#plot bands and DOS with no SOC
calcs_tb.plot.opdos(en_range=[-4,4])
calcs_tb.plot.apdos(en_range=[-4,4])
calcs_tb.plot.bands(DOSPlot="APDOS",en_range=[-4,4])

# change to fully relativistic pseudos
calcs.change_pseudos("./FR_PP")

# add SOC keywords to input files
calcs.change_input("&system","lspinorb",".true.")
calcs.change_input("&system","noncolin",".true.")

# remove nspin variable from inputs
calcs.change_input("&system","nspin",None)

# generate PAO-TB hamiltonian with SOC
calcs_tb=calcs.tight_binding(proj_thr=0.95,tb_kp_factor=4,
                             symmetrize=True,sym_thr=1.e-8,
                             sym_max_iter=50)

#calculate bands,DOS, AHC, and SHC
calcs_tb.dos()
calcs_tb.bands()
calcs_tb.ahc()
calcs_tb.shc()

#plot bands and DOS with SOC
calcs_tb.plot.opdos(en_range=[-4,4])
calcs_tb.plot.apdos(en_range=[-4,4])
calcs_tb.plot.bands(DOSPlot="APDOS",en_range=[-4,4])

# plot AHC,SHC, mag dichro, spin dichro
calcs_tb.plot.shc(en_range=[-4,4])
calcs_tb.plot.ahc(en_range=[-4,4])

# run/submit the calcs
calcs.submit()



