import sys, subprocess
import PySide
import matplotlib
matplotlib.use('Agg')
import fitsubs
import numpy as np


def usage():
	print "Gaussian fitting of Pseudo-Atomic Orbitals (PAOs)"
	S = '''
		Usage: gauss_pao.py <atom label> <path of PP>

		atom label 	- atom label
		path of PP 	- path of Norm-conserving PP
	'''
	print S

def get_atno(atomlabel) :

	mass = {
        
	"H": (1, 1.008, 'Hydrogen'),
	"He": (2, 4.003, 'Helium'),
	"Li": (3, 6.938, 'Lithium'),
	"Be": (4, 9.012, 'Beryllium'),
	"B": (5, 10.806, 'Boron'),
	"C": (6, 12.010, 'Carbon'),
	"N": (7, 14.006, 'Nitrogen'),
	"O": (8, 15.999, 'Oxygen'),
	"F": (9, 18.998, 'Fluorine'),
	"Ne": (10, 20.180, 'Neon'),
	"Na": (11, 22.990, 'Sodium'),
	"Mg": (12, 24.304, 'Magnesium'),
	"Al": (13, 26.982, 'Aluminium'),
	"Si": (14, 28.084, 'Silicon'),
	"P": (15, 30.974, 'Phosphorus'),
	"S": (16, 32.059, 'Sulfur'),
	"Cl": (17, 35.446, 'Chlorine'),
	"Ar": (18, 39.948, 'Argon'),
	"K": (19, 39.098, 'Potassium'),
	"Ca": (20, 40.078, 'Calcium'),
	"Sc": (21, 44.955, 'Scandium'),
	"Ti": (22, 47.867, 'Titanium'),
	"V": (23, 50.941, 'Vanadium'),
	"Cr": (24, 51.996, 'Chromium'),
	"Mn": (25, 54.938, 'Manganese'),
	"Fe": (26, 55.845, 'Iron'),
	"Co": (27, 58.933, 'Cobalt'),
	"Ni": (28, 58.693, 'Nickel'),
	"Cu": (29, 63.546, 'Copper'),
	"Zn": (30, 65.380, 'Zinc'),
	"Ga": (31, 69.723, 'Gallium'),
	"Ge": (32, 72.630, 'Germanium'),
	"As": (33, 74.922, 'Arsenic'),
	"Se": (34, 78.971, 'Selenium'),
	"Br": (35, 79.901, 'Bromine'),
	"Kr": (36, 83.798, 'Krypton'),
	"Rb": (37, 85.468, 'Rubidium'),
	"Sr": (38, 87.620, 'Strontium'),
	"Y": (39, 88.906, 'Yttrium'),
	"Zr": (40, 91.224, 'Zirconium'),
	"Nb": (41, 92.906, 'Niobium'),
	"Mo": (42, 95.950, 'Molybdenum'),
    	"Tc": (43, 98.0, 'Tecnetium'),
	"Ru": (44, 101.070,' Ruthenium'),
	"Rh": (45, 102.906, 'Rhodium'),
	"Pd": (46, 106.420, 'Palladium'),
	"Ag": (47, 107.868, 'Silver'),
	"Cd": (48, 112.414, 'Cadmium'),
	"In": (49, 114.818, 'Indium'),
	"Sn": (50, 118.710, 'Tin'),
	"Sb": (51, 121.760, 'Antimony'),
	"Te": (52, 127.600, 'Tellurium'),
	"I": (53, 126.904, 'Iodine'),
	"Xe": (54, 131.293, 'Xenon'),
	"Cs": (55, 132.905, 'Cesium'),
	"Ba": (56, 137.327, 'Barium'),
	"La": (57, 138.905, 'Lanthanum'),
	"Ce": (58, 140.116, 'Cerium'),
	"Pr": (59, 140.908, 'Praseodymium'),
	"Nd": (60, 144.242, 'Neodymium'),
	"Sm": (62, 150.360, 'Samarium'),
	"Eu": (63, 151.964, 'Europium'),
	"Gd": (64, 157.250, 'Gadolinium'),
	"Tb": (65, 158.925, 'Terbium'),
	"Dy": (66, 162.500, 'Dysprosium'),
	"Ho": (67, 164.930, 'Holmium'),
	"Er": (68, 167.259, 'Erbium'),
	"Tm": (69, 168.934, 'Thulium'),
	"Yb": (70, 173.054, 'Ytterbium'),
	"Lu": (71, 174.967, 'Lutetium'),
	"Hf": (72, 178.490, 'Hafnium'),
	"Ta": (73, 180.947, 'Tantalum'),
	"W": (74, 183.840, 'Tungsten'),
	"Re": (75, 186.207, 'Rhenium'),
	"Os": (76, 190.230, 'Osmium'),
	"Ir": (77, 192.217, 'Iridium'),
	"Pt": (78, 195.084, 'Platinum'),
	"Au": (79, 196.966, 'Gold'),
	"Hg": (80, 200.592, 'Mercury'),
	"Tl": (81, 204.382, 'Thallium'),
	"Pb": (82, 207.200, 'Lead'),
	"Bi": (83, 208.980, 'Bismuth'),
    	"Po": (84, 209, 'Polonium'),
	"At": (85, 210, 'Astatine'),
	"Th": (90, 232.038, 'Thorium'),
	"Pa": (91, 231.036, 'Protactinium'),
	"U": (92, 238.029, 'Uranium'),
	"Am": (95, 243.000, 'Americium'),}

	try:
		atno = mass[atomlabel][0]
		return atno
	except Exception, e:
		print "Atom is not in List provide atom. no."
		atno = input("atom no.?")
		return int(float(atno))

	
def fitPAOs(atomlabel,pp_path):
	
	import matplotlib.pyplot as plt
	import os,re

	gbs_path = "%s_sto3g.gbs"%atomlabel
	if not os.path.exists(gbs_path):
		print 'STO-3G file not foundi - Please include file in current working directory'
		sys.exit()

	str="grep CHI %s|grep label|awk '{print $6}'|cut -d'\"' -f2|paste -s"%pp_path
	orb_lbls=subprocess.check_output(str,shell=True).split()
	print orb_lbls
	gbs_lines = open(gbs_path,'r').read()
	r1 = re.compile(r'[a-zA-Z]+\s*\d+',re.MULTILINE)
	block_lbls = r1.findall(gbs_lines)

	for i in range(len(orb_lbls)):
		ll = i
		orb_lbl = orb_lbls[i][-1]
		if orb_lbl.lower() == 's':l = 0
		elif orb_lbl.lower() == 'p':l = 1
		elif orb_lbl.lower() == 'd':l = 2

		for j in range(len(block_lbls)):
			if orb_lbl.lower() in block_lbls[j] or orb_lbl.upper() in block_lbls[j]:
				nblock = j+1
				print orb_lbl, block_lbls[j],nblock
				break

		nent = 3 #No. of entries in block - same for all sto3g files - 3 Gaussians

		subblocks = [nblock]


		#Loads the initial guess Gaussians.
		prim_gauss = fitsubs.get_prim_gauss(atomlabel,gbs_path,subblocks)
		#Loads the wfcs from the UPF file
		rmesh, wfc_mat = fitsubs.read_qe_radial_function(pp_path)
		l_phi = wfc_mat[ll]

		#the standard convention s,p,d,f = 0,1,2,3 
		l_index = prim_gauss[:,0] == l
		prim_gauss_cart = prim_gauss[l_index,:] #these are not normalized




		#Fitting methods
		#==============Method-1:optimized gaussians: coefficients=======================================#
		opt1_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph_coeffs_only(prim_gauss_cart,l_phi,rmesh,l)
		#==============Method-2:optimized gaussians: exponents and coefficients from opt coeffs=========#
		opt2_prim_gauss_sph = np.concatenate(( np.zeros((nent,1)),opt1_prim_gauss_sph ),axis=1)
		opt3_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph(opt2_prim_gauss_sph,l_phi,rmesh,l)
		#==============Method-3:optimized gaussians: exponents and coefficientsi========================#
		opt_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph(prim_gauss_cart,l_phi,rmesh,l)

		print "##################################################"
		print "prim_gauss_cart    :\n"  , prim_gauss_cart 
		print "opt1_prim_gauss_sph:\n"  , opt1_prim_gauss_sph
		print "opt3_prim_gauss_sph:\n"  , opt3_prim_gauss_sph
		print "opt_prim_gauss_sph :\n"  , opt_prim_gauss_sph


		os.system('clear')
		print "==============cartesian gaussians"
		title='cartesian gaussians'
		fitsubs.plot_prim_gauss(prim_gauss_cart,rmesh,l,l_phi,title)
		plt.ylim((-0.5, 1.5))
		plt.draw()
		plt.savefig("%s_%s_fig1.pdf"%(atomlabel,orb_lbl.upper()))

		print "==============optimized gaussians: coefficients"
		title='opt gauss: coefficients only'
		fitsubs.plot_prim_gauss_sph(opt1_prim_gauss_sph,rmesh,l,l_phi,title)
		plt.ylim((-0.5, 1.5))
		plt.draw()
		plt.savefig("%s_%s_fig2.pdf"%(atomlabel,orb_lbl.upper()))

		print "==============optimized gaussians: exponents and coefficients from opt coeffs"
		title='opt gauss: exponents and coefficients from opt coeffs'
		print title
		fitsubs.plot_prim_gauss_sph(opt3_prim_gauss_sph,rmesh,l,l_phi,title)
		plt.ylim((-0.5, 1.5))
		plt.draw()
		plt.savefig("%s_%s_fig3.pdf"%(atomlabel,orb_lbl.upper()))

		print "==============optimized gaussians: exponents and coefficients"
		title='opt gauss: exponents and coefficients'
		fitsubs.plot_prim_gauss_sph(opt_prim_gauss_sph,rmesh,l,l_phi,title)
		plt.ylim((-0.5, 1.5))
		plt.draw()
		plt.savefig("%s_%s_fig4.pdf"%(atomlabel,orb_lbl.upper()))


		fullpath='./%s_atom'%atomlabel
		if not os.path.exists(fullpath):
			os.system('mkdir %s'%fullpath)

		##==================UNCOMMENT TO CHOOSE FITTING METHOD=================================
		#Generatig fit with Method-3
		fitsubs.print_opt_prim_gauss_cart(opt_prim_gauss_sph,fullpath,atomlabel,l)
		#Generatig fit with Method-2
		#fitsubs.print_opt_prim_gauss_cart(opt3_prim_gauss_sph,fullpath,atomlabel,l)

def mk_basis(atomlabel,pp_path):

	
	import sys,commands
	from string import digits
	import fitsubs

	#basis_file 
	qespresso_order={0:[0],1:[0,1,-1],2:[0,1,-1,2,-2]}
	#follows UPF order
	atno      = get_atno(atomlabel)

	ppNm = pp_path

	orbs = commands.getoutput("grep CHI %s |grep label|cut -d= -f6|awk '{print $1}'"%ppNm).translate(None,digits).replace("\"","").split()
	basis_file = []
	basis_type = []
	for i in orbs:
		basis_file.append(atomlabel+"."+i.lower()+".py")
		if i == "s" or i == "S": basis_type.append(0)
		elif i == "p" or i == "P": basis_type.append(1)
		elif i == "d" or i == "D": basis_type.append(2)
		else: print "Unidentified orbital"

	basis_file_path = './%s_atom'%atomlabel

	f = open(basis_file_path+'/'+'%s_basis.py'%atomlabel,'w+')
	f.write('basis_data={ %d :[\n' % atno)
	total_f = len(basis_file)
	print 'Total number of files: %d\n' % total_f
	for fcounter, ibf in enumerate(basis_file):
	    print 'Reading:'+basis_file_path+'/'+ibf
	    basisFile = basis_file_path+'/'+ibf
	    execfile(basis_file_path+'/'+ibf,globals())
	    il = basis_type[fcounter]
	    total_m = len(qespresso_order[il])
	    for mcounter,im in enumerate(qespresso_order[il]):
		npgto = len( cgto[(basis_type[fcounter],im)] )
		for counter,ipgto in enumerate(cgto[(basis_type[fcounter],im)]):
		    if counter == 0:
		       if mcounter == 0:
			  ssymb = '[['
		       else:
			  ssymb = ' ['
		       esymb = ',\n'
		    elif counter == npgto-1:
		       ssymb = '  '
		       if mcounter == total_m-1:
			  if fcounter == total_f-1:
			     esymb = ']]\n'
			  else:
			     esymb = ']],\n'
		       else:
			  esymb = '],\n'
		    else:
		       ssymb = '  '
		       esymb = ',\n'
		    f.write('%s(%d, %d, %d, %17.10f, %17.10f)%s'%(ssymb,ipgto[0],ipgto[1],ipgto[2],ipgto[3],ipgto[4],esymb))
	f.write(']}\n')
	f.close()

if __name__ == '__main__':

	if len(sys.argv) > 2:

		atomlabel=sys.argv[1]
		pp_path = sys.argv[2]
		#fitting 
		fitPAOs(atomlabel,pp_path)
		mk_basis(atomlabel,pp_path)
	else:
		print 'Error in No of Arguments'
		usage()
		sys.exit()

