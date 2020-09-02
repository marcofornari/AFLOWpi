#!/usr/bin/env python
#######################################################################
# Fit UPF radial pseudowavefunctions with gaussian orbitals
# Davide Ceresoli - May 2016
#
# Notes:
# - UPFv1 files must be embedded in <UPF version="1.0">...</UPF> element
# - contraction coefficients for d and f orbitals correspond to the
#   cubic harmonics
#######################################################################
import sys
import numpy as np
import argparse
import io
from math import *
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from xml.etree import ElementTree as ET
import  matplotlib
matplotlib.use("pdf")
from matplotlib import pylab
import os

def get_atom_no(atom_n):
    spn_map={"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,"Na":11,
             "Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,"K":19,"Ca":20,
             "Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,"Cu":29,
             "Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,"Rb":37,
             "Sr":38,"Y":39,"Zr":40,"Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,
             "Ag":47,"Cd":48,"In":49,"Sn":50,"Sb":51,"Te":52,"I":53,"Xe":54,"Cs":55,
             "Ba":56,"La":57,"Ce":58,"Pr":59,"Nd":60,"Pm":61,"Sm":62,"Eu":63,
             "Gd":64,"Tb":65,"Dy":66,"Ho":67,"Er":68,"Tm":69,"Yb":70,"Lu":71,
             "Hf":72,"Ta":73,"W":74,"Re":75,"Os":76,"Ir":77,"Pt":78,"Au":79,"Hg":80,
             "Tl":81,"Pb":82,"Bi":83,"Po":84,"At":85,"Rn":86,"Fr":87,"Ra":88,
             "Ac":89,"Th":90,"Pa":91,"U":92,"Np":93,"Pu":94,"Am":95,"Cm":96,"Bk":97,
             "Cf":98,"Es":99,"Fm":100,"Md":101,"No":102,"Lr":103,"Rf":104,"Db":105,
             "Sg":106,"Bh":107,"Hs":108,"Mt":109,"Ds":110,"Rg":111,"Cn":112}

    return spn_map[atom_n]
    
    

# double factorial (n!!)
def fact2(n):
    if n <= 1: return 1
    return n*fact2(n-2)


#======================================================================
# GTO orbital
#======================================================================
def gto(r, l, params):
    alpha, beta = params[0:2]
    coeffs = params[2:]

    gto = np.zeros_like(r)
    for (j,coeff) in enumerate(coeffs):
        #assert alpha > 0 and beta > 0
        zeta = alpha/beta**j
        i = np.where(zeta*r*r > -12.0)
        gto[i] += r[i]**l * coeff*np.exp(-zeta*r[i]*r[i])

    return gto


#======================================================================
# Target function whose least square has to be minimized
#======================================================================
def target(params, r, rab, wfc, l):
    return wfc - r*gto(r, l, params)

def target_squared(params, r, rab, wfc, l):
    diff = target(params, r, rab, wfc, l)
    return np.dot(diff, diff)

#======================================================================
# Fit radial wfc with gaussians
#======================================================================
def fit(nzeta, label, l, r, rab, wfc,threshold):
    assert len(wfc) == len(r)

    wfc = np.array(wfc)
    r = np.array(r)

    params0 = np.array([4.0, 4.0]) # initial alpha and beta
    params0 = np.append(params0, np.ones((nzeta,)))

    # least squares
    if True:
        params, fit_cov, info, mesg, ier = \
            leastsq(target, params0, args=(r, rab, wfc, l), full_output=1, \
            maxfev=50000, ftol=1e-10, xtol=1e-10)
        if ier > 0:
            print(("ERROR: ier=", ier, "mesg=", mesg))
            print(("ERROR: info[nfev]=", info["nfev"]))
            print(("ERROR: info[fvec]=", sum(info["fvec"]**2.0)))
    else: # minimize
        opt = minimize(target_squared, params0, args=(r, rab, wfc, l), \
                       method='CG', tol=1e-10)
        #opt = basinhopping(target_squared, params0, minimizer_kwargs={'args':(r, rab, wfc, l)})
        params = opt.x
        if not opt.success:
           print(("ERROR: opt.status=", opt.status))
           print(("ERROR: opt.message=", opt.message))
           print(("ERROR: opt.nfev=", opt.nfev))
           print(("ERROR: opt.fun=", opt.fun))


    alpha, beta = params[0:2]
    n = sqrt(fact2(2*l+1)/(4.0*pi))
    coeffs = params[2:] * n
    expon = []
    print("alpha = %f, beta = %f" % (alpha, beta))
    for (j,coeff) in enumerate(coeffs):
        zeta = alpha/beta**j
        expon.append(zeta)
        print("coeff = %f,  zeta = %f" % (coeff, zeta))

    with open("wfc"+label+".dat", "w") as f:
        gto_r = gto(r, l, params)
    #    for i in range(len(wfc)):
    #        f.write("%f %f %f\n" % (r[i], wfc[i], r[i]*gto_r[i]))
    pylab.plot(r, wfc, '.', label=label+"_orig")
    pylab.plot(r, r*gto(r, l, params), label=label+"_fit")
    pylab.xlim(0,r[-1])
    res=target_squared(params, r, rab, wfc, l)
    print("INFO: fit result:", res)

    quit_flag=False
    if np.abs(res)>threshold:
#        print("exiting")
        quit_flag=True


    
    return coeffs, expon,quit_flag



#======================================================================
# Print python block for orbitals
#======================================================================
def print_python_block(bfile, label, l, coeffs, expon):
    nzeta = len(coeffs)

    ret_str=''
    ret_str+= "# label= %s l= %s\n"%(label,l)
    if l == 0:
        ret_str+= "[[\n"
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,0,0,coeffs[n],expon[n])
        ret_str+= "]],\n"

    elif l == 1:
        ret_str+= "[[\n"
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,0,1,coeffs[n],expon[n])
        ret_str+= "], [\n"
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,1,0,coeffs[n],expon[n])
        ret_str+= "], [\n"
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,0,0,coeffs[n],expon[n])
        ret_str+= "]],\n"

    elif l == 2:
        ret_str+= "[[\n"

        fact = 0.5/sqrt(3.0)
        for n in range(nzeta):  # 1/(2.0*sqrt(3))*(2*z2 - x2 - y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,0,2,2.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,2,0,-fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (2,0,0,-fact*coeffs[n],expon[n])
        ret_str+= "], [\n"

        for n in range(nzeta): # xz
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,0,1,coeffs[n],expon[n])
        ret_str+= "], [\n"

        for n in range(nzeta): # yz
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,1,1,coeffs[n],expon[n])
        ret_str+= "], [\n"

        fact = 0.5
        for n in range(nzeta): # 1/2 * (x2 - y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (2,0,0,fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,2,0,-fact*coeffs[n],expon[n])
        ret_str+= "], [\n"

        for n in range(nzeta): # xy
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,1,0,coeffs[n],expon[n])
        ret_str+= "]],\n"

    elif l == 3:        
        # fz3,fxz2,fyz2,fz(x2-y2),fxyz,fx(x3-3y2),fy(3x2-y2)

        ret_str+= "[[\n"

        fact=0.5/np.sqrt(15)
        for n in range(nzeta):  # 1/(2*sqrt(15)) * z*(2*z2 - 3*x2 - 3*y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,0,3, 2.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (2,0,1,-3.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,2,1,-3.0*fact*coeffs[n],expon[n])
        ret_str+= "], [\n"
        fact=0.5/np.sqrt(10)
        for n in range(nzeta):  # 1/(2*sqrt(10)) * x*(4*z2 - x2 - y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,0,2, 4.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,0,3,-1.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,2,0,-1.0*fact*coeffs[n],expon[n])
        ret_str+= "], [\n"
        fact=0.5/np.sqrt(10)
        for n in range(nzeta):  # 1/(2*sqrt(10)) * y*(4*z2 - x2 - y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,1,2, 4.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (2,1,0,-1.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,3,0,-1.0*fact*coeffs[n],expon[n])
        ret_str+= "], [\n"
        fact=0.5           
        for n in range(nzeta):  # 1/2 * z*(x2 - y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (2,0,1,     fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,2,1,-1.0*fact*coeffs[n],expon[n])
        ret_str+= "], [\n"
        fact=1.0
        for n in range(nzeta):  # x*y*z
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,1,1,fact*coeffs[n],expon[n])
        ret_str+= "], [\n"
        fact=0.5/np.sqrt(6) 
        for n in range(nzeta):  # 1/(2*sqrt(6) * x*(x2 - 3*y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (3,0,0,     fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (1,2,0,-3.0*fact*coeffs[n],expon[n])
        ret_str+= "], [\n"
        fact=0.5/np.sqrt(6) 
        for n in range(nzeta):  # 1/(2*sqrt(6) * y*(3*x2 - y2)
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (2,1,0, 3.0*fact*coeffs[n],expon[n])
        for n in range(nzeta):
            ret_str+= "   (%i,%i,%i,%20.10f,%20.10f),\n" % (0,3,0,-1.0*fact*coeffs[n],expon[n])
        ret_str+= "]],\n"

    return ret_str




def gen_gauss_proj(nzeta,xml_file,threshold=0.5,exclude=''):

    print(xml_file)
    base_dir=os.path.dirname(xml_file)
    #### open the UPF file ####
    try:
        with open(xml_file) as f:
            xml_file_content = f.read()
        xml_file_content = xml_file_content.replace('&', '&amp;')
        root = ET.fromstring(xml_file_content)
        version = root.attrib["version"]
        upfver = int(version.split(".")[0])
        print("INFO: fitting file", xml_file, "with", nzeta, "gaussians")
        print("INFO: UPF version", upfver, "detected")
    except Exception as inst:
        try:
            version = root.attrib["UPF version"]
            upfver = int(version.split(".")[0])
            print("INFO: fitting file", xml_file, "with", nzeta, "gaussians")
            print("INFO: UPF version", upfver, "detected")
        except Exception as inst:
            print("Unexpected error opening %s: %s" % (xml_file, inst))
            sys.exit(1)


    #### get element name ####
    if upfver == 1:
        text = root.find('PP_HEADER').text.split()
        for i in range(len(text)):
            if text[i] == 'Element':
                element = text[i-1].strip()
                break
    else:
        element = root.find('PP_HEADER').attrib["element"].strip()

    atno = get_atom_no(element)    

    print("INFO: element=", element, "atomic number=", atno)
    print()


    #### open basis file ####

    #basisfile.write()
    pylab.title(xml_file)


    #### read the radial grid ####
    text = root.find('PP_MESH/PP_R').text
    r = np.array( [float(x) for x in text.split()] )
    text = root.find('PP_MESH/PP_RAB').text
    rab = np.array( [float(x) for x in text.split()] )

    basisfile=''
    #### read and fit radial wavefunctions ####
    if upfver == 1:
        pot = root.find('PP_LOCAL')
        if pot is None: quit()
        v = [float(x) for x in pot.text.split()]
        #f = open('vlocal.dat', 'w')
        #for i in range(len(v)):
        #    f.write("%f %f\n" % (r[i], v[i]))
        #f.close()

        chis = root.find('PP_PSWFC')
        if chis is None:
             print("ERROR: cannot find PP_PSWFC tag")
             sys.exit(1)
        data = io.StringIO(chis.text)
        nlines = len(r)/4
        if len(r) % 4 != 0: nlines += 1

        while True:
            line = data.readline()
            if line == "\n": continue
            if line == "": break
            label, l, occ, dummy = line.split()
            l = int(l)
            occ = float(occ)
            wfc = []

            for i in range(nlines):
                wfc.extend(list(map(float, data.readline().split())))
            wfc = np.array(wfc)

            if exclude.find(label) >= 0:
                print("INFO: skipping", label)
                continue

            norm = sum(wfc*wfc*rab)
            print("INFO: fitting pswfc", label, "l=", l, "norm=", norm)
            #wfc *= 1.0/sqrt(norm)
            coeffs, expon,quit_flag = fit(nzeta, label, l, r, rab, wfc,threshold)
            ret_str += print_python_block(basisfile, label, l, coeffs, expon)
            print()

            if quit_flag:
                raise SystemExit

        basisfile = open(os.path.join(base_dir,element+'_basis.py'), "w")        
        ret_str="basis_data = { %i : [\n" % (atno) + ret_str

        basisfile.write(ret_str)

        betas = root.find('PP_NONLOCAL/PP_BETA')
        if betas is None:
             print("ERROR: cannot find PP_BETA tag")
             sys.exit(1)
        data = io.StringIO(betas.text)

        while True:
            line = data.readline()
            if line == "\n": continue
            if line.strip() == "": break
            ibeta, l, dummy, dummy = line.split()
            ibeta = int(ibeta)
            l = int(l)

            npoints = int(data.readline())
            nlines = npoints/4
            if npoints % 4 != 0: nlines += 1
            beta = []

            for i in range(nlines):
                beta.extend(list(map(float, data.readline().split())))
            print(beta)
            beta = np.array(beta)
            #f = open("beta_%i_%i.dat" % (ibeta, l), 'w')
            #for i in range(len(beta)):
            #    f.write("%f %f\n" % (r[i], beta[i]))
            #f.close()
            line = data.readline()
            #print(line)
            line = data.readline()
            #print(line)

    else:
        pot = root.find('PP_LOCAL')
        if pot is None: quit()
        v = [float(x) for x in pot.text.split()]
        #f = open('vlocal.dat', 'w')
        #for i in range(len(v)):
        #    f.write("%f %f\n" % (r[i], v[i]))
        #f.close()

        #pot = root.find('PP_GIPAW/PP_GIPAW_VLOCAL/PP_GIPAW_VLOCAL_AE')
        #if pot is None: quit()
        #v = [float(x) for x in pot.text.split()]
        #f = open('vlocal_ae.dat', 'w')
        #for i in xrange(len(v)):
        #    f.write("%f %f\n" % (r[i], v[i]))
        #f.close()

        i = 0

        ret_str=''
        while True:
            i += 1
            chi = root.find('PP_PSWFC/PP_CHI.%i' % (i))
            if chi is None: break

            label = chi.attrib["label"]
            if exclude.find(label) >= 0:
                print("INFO: skipping", label)
                continue
            l = int(chi.attrib["l"])
            wfc = [float(x) for x in chi.text.split()]
            assert len(wfc) == len(r)
            wfc = np.array(wfc)

            norm = sum(wfc*wfc*rab)
            print("INFO: fitting pswfc", label, "l=", l, "norm=", norm)
            #wfc *= 1.0/sqrt(norm)
            coeffs, expon,quit_flag = fit(nzeta, label, l, r, rab, wfc,threshold)

            ret_str += print_python_block(basisfile, label, l, coeffs, expon)

            if quit_flag:
                return quit_flag


        basisfile = open(os.path.join(base_dir,element+'_basis.py'), "wt")        
        ret_str="basis_data = { %i : [\n" % (atno) + ret_str

        basisfile.write(ret_str)        
    #    basisfile.write(ret_str)

    basisfile.write("]}\n")
    basisfile.close()
    # pylab.legend(loc='lower right')
    # pylab.grid()
    # pylab.xlabel('r (bohrradius)')
    # pylab.ylabel('radial wfc')
    # pylab.show()

    # pfn= basisfile.name.split('_')[0].split('.')[0].split('-')[0]+'.pdf'
    # pylab.savefig(pfn)
    print("INFO: file", basisfile.name, "created!")

    return quit_flag



def gauss_fit_loop(xml_file,threshold=0.5):
    nzeta=2

    while True:
        qf = gen_gauss_proj(nzeta,xml_file,exclude='',threshold=threshold)
        if qf==False:
            break
        else:
            nzeta+=1
    





