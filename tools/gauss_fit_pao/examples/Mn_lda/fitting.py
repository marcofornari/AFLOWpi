#By Luis Agapito at Marco Buongiorno Nardelli's UNT group
#2013

import matplotlib.pyplot as plt
import sys
from sympy import *
from __future__ import division
#

#tt = t tilde
# table of coefficients
Nnm = {      
(0, 0) : 1/4,
(1, 0) : 1/4,
(1, 1) : 1/4,
(2, 0) : 3  ,
(2, 1) : 1/4,
(2, 2) : 1  ,
(2,-2) : 1/4,
(3, 0) : 15 ,
(3, 1) : 10 ,
(3, 2) : 1  ,
(3,-2) : 1/4,
(3, 3) : 6   
}
#Conversion Table. Appendix B of Mathar's paper
def f000():
 return ttt(0,0,0)

def f001():
 return ttt(1,0,0)

def f100():
 return ttt(1,1,0)

#added
def f010():
 return ttt(1,-1,0)

def f002():
 return 1/3*(ttt(2,0,0)+ttt(0,0,2))

def f020():
 return -1/6*ttt(2,0,0) - 1/2*ttt(2,2,0) + 1/3*ttt(0,0,2)

#added
def f200():
 return -1/6*ttt(2,0,0) + 1/2*ttt(2,2,0) + 1/3*ttt(0,0,2)

def f011():
 return ttt(2,-1,0)

#added
def f101():
 return ttt(2,1,0)

def f110():
 return ttt(2,-2,0)

def f003():
 return 1/5*ttt(3,0,0) + 3/5*ttt(1,0,2)

def f030():
 return -1/4*ttt(3,-3,0) - 3/20*ttt(3,-1,0) + 3/5*ttt(1,-1,2)

#added
def f300():
 return  1/4*ttt(3,3,0) - 3/20*ttt(3,1,0) + 3/5*ttt(1,1,2) 

def f012():
 return 1/5*(ttt(3,-1,0) + ttt(1,-1,2))

#added
def f102():
 return 1/5*( ttt(3,1,0) + ttt(1,1,2) )

def f021():
 return -1/10*ttt(3,0,0) - 1/2*ttt(3,2,0) + 1/5*ttt(1,0,2)

#added
def f201():
 return -1/10*ttt(3,0,0) + 1/2*ttt(3,2,0) + 1/5*ttt(1,0,2)

def f120():
 return -1/4*ttt(3,3,0) - 1/20*ttt(3,1,0) + 1/5*ttt(1,1,2)

#added
def f210():
 return  1/4*ttt(3,-3,0) -1/20*ttt(3,-1,0) + 1/5*ttt(1,-1,2) 

def f111():
 return ttt(3,-2,0)

#cartesian2spherical                 
c2s= {
(0,0,0): f000,
(0,0,1): f001,
(1,0,0): f100,
(0,1,0): f010,
(0,0,2): f002,
(0,2,0): f020,
(2,0,0): f200,
(0,1,1): f011,
(1,0,1): f101,
(1,1,0): f110,
(0,0,3): f003,
(0,3,0): f030,
(3,0,0): f300,
(0,1,2): f012,
(1,0,2): f102,
(0,2,1): f021,
(2,0,1): f201,
(1,2,0): f120,
(2,1,0): f210,
(1,1,1): f111
}

s  = Symbol('s')
pz = Symbol('pz') 
px = Symbol('px') 
py = Symbol('py') 
dz2 = Symbol('dz2') 
dzx = Symbol('dzx') 
dzy = Symbol('dzy') 
dx2_y2 = Symbol('dx2_y2') 
dxy = Symbol('dxy') 

r  = Symbol('r')

# From QE projwfc.x definitions
# for l=1:
#   1 pz     (m=0)
#   2 px     (real combination of m=+/-1 with cosine)
#   3 py     (real combination of m=+/-1 with sine)
# 
# for l=2:
#   1 dz2    (m=0)
#   2 dzx    (real combination of m=+/-1 with cosine)
#   3 dzy    (real combination of m=+/-1 with sine)
#   4 dx2-y2 (real combination of m=+/-2 with cosine)
#   5 dxy    (real combination of m=+/-2 with sine)

# From http://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics_with_l_.3D_4
#            Y^{m,c/s}_n
#   1 pz     Y^{0,c}_1
#   2 px     Y^{1,c}_1
#   3 py     Y^{1,s}_1
#
#   1 dz2    Y^{0,c}_2
#   2 dzx    Y^{1,c}_2
#   3 dzy    Y^{1,s}_2
#   4 dx2-y2 Y^{2,c}_2
#   5 dxy    Y^{2,s}_2

#( n, m, 'c/s')
basis_name = {
 ( 0, 0, 'c') : s  ,
 ( 1, 0, 'c') : pz ,
 ( 1, 1, 'c') : px ,
 ( 1, 1, 's') : py ,
 ( 2, 0, 'c') : dz2 ,
 ( 2, 1, 'c') : dzx ,
 ( 2, 1, 's') : dzy ,
 ( 2, 2, 'c') : dx2_y2 ,
 ( 2, 2, 's') : dxy 
}

 
#if not (lx,ly,lz) in mydict
# sys.exit("The definition of a cartesian Gaussian with lx,ly,lz is not defined. Exiting ...")
#
#f_expan = mydict[(lx,ly,lz)]()

def Mathar_Y(n,m):
 global basis_name
 if m >= 0:
  symbol = basis_name[(n,m,'c')]
 else:
  symbol = basis_name[(n,-m,'s')]
 #print symbol
 return symbol

def get_Nnm(n,m):
 global Nnm
 if (m%2 != 0) & (m<0) :
  return Nnm[(n,-m)]
 else:
  return Nnm[(n,m)]

def tt(n,m):
 import numpy as np
 import scipy.misc as spm
 tt_expression=4*sqrt(np.pi*get_Nnm(n,m)/spm.factorial2(2*n+1))*r**n*exp(-1*r**2)*Mathar_Y(n,m)
 return tt_expression

#3-index definition of tt. Appendix B  of Mathar's paper
def ttt(n,m,s):
 return (r**s)*tt(n,m)

##################################################
#tests
c2s[(0,0,0)]()

sqrt(3/(4*np.pi))*c2s[(1,0,0)]()
sqrt(3/(4*np.pi))*c2s[(0,1,0)]()
sqrt(3/(4*np.pi))*c2s[(0,0,1)]()

1/4*sqrt(5/np.pi)*simplify(-c2s[(2,0,0)]()- c2s[(0,2,0)]() + 2*c2s[(0,0,2)]() )
1/2*sqrt(15/np.pi)*simplify(c2s[(0,1,1)]())
1/2*sqrt(15/np.pi)*simplify(c2s[(1,0,1)]())
1/2*sqrt(15/np.pi)*simplify(c2s[(1,1,0)]())
1/4*sqrt(15/np.pi)*simplify(c2s[(2,0,0)]()- c2s[(0,2,0)]()) 


c2s[(0,0,0)]()
c2s[(1,0,0)]()
c2s[(0,1,0)]()
c2s[(0,0,1)]()
c2s[(1,1,0)]()
c2s[(1,0,1)]()
c2s[(0,1,1)]()

#this results are correct
#>>> c2s[(0,0,0)]()
#3.54490770181103*s*exp(-r**2)
#>>> 
#>>> sqrt(3/(4*np.pi))*c2s[(1,0,0)]()
#1.0*px*r*exp(-r**2)
#>>> sqrt(3/(4*np.pi))*c2s[(0,1,0)]()
#1.0*py*r*exp(-r**2)
#>>> sqrt(3/(4*np.pi))*c2s[(0,0,1)]()
#1.0*pz*r*exp(-r**2)
#>>> 
#>>> 1/4*sqrt(5/np.pi)*simplify(-c2s[(2,0,0)]()- c2s[(0,2,0)]() + 2*c2s[(0,0,2)]() )
#1.0*dz2*r**2*exp(-r**2)
#>>> 1/2*sqrt(15/np.pi)*simplify(c2s[(0,1,1)]())
#1.0*dzy*r**2*exp(-r**2)
#>>> 1/2*sqrt(15/np.pi)*simplify(c2s[(1,0,1)]())
#1.0*dzx*r**2*exp(-r**2)
#>>> 1/2*sqrt(15/np.pi)*simplify(c2s[(1,1,0)]())
#1.0*dxy*r**2*exp(-r**2)
#>>> 1/4*sqrt(15/np.pi)*simplify(c2s[(2,0,0)]()- c2s[(0,2,0)]()) 
#1.0*dx2_y2*r**2*exp(-r**2)

#>>> c2s[(0,0,0)]()
#3.54490770181103*s*exp(-r**2)
#>>> c2s[(1,0,0)]()
#2.04665341589298*px*r*exp(-r**2)
#>>> c2s[(0,1,0)]()
#2.04665341589298*py*r*exp(-r**2)
#>>> c2s[(0,0,1)]()
#2.04665341589298*pz*r*exp(-r**2)
#>>> c2s[(1,1,0)]()
#0.915291232863769*dxy*r**2*exp(-r**2)
#>>> c2s[(1,0,1)]()
#0.915291232863769*dzx*r**2*exp(-r**2)
#>>> c2s[(0,1,1)]()
#0.915291232863769*dzy*r**2*exp(-r**2)
#>>> 


##################################################
#reading a Gaussian-format basis set

import fitsubs

atomlabel='Si'
gbs_path ='/Users/believe/Google\ Drive/unt2/xintegrals/basis/631.gbs'
subblocks = [3,4]

prim_gauss = get_prim_gauss(atomlabel,gbs_path,subblocks)

atomlabel='H'
gbs_path ='/Users/believe/Google\ Drive/unt2/xintegrals/basis/321.gbs'
subblocks = [1,2]
prim_gauss = get_prim_gauss(atomlabel,gbs_path,subblocks)
plot_contracted_gauss(prim_gauss)

atomlabel='H'
gbs_path ='/Users/believe/Google\ Drive/unt2/xintegrals/basis/sto3g.gbs'
subblocks = [1]
prim_gauss = get_prim_gauss(atomlabel,gbs_path,subblocks)
#plot_contracted_gauss(prim_gauss,0)
plot_contracted_gauss_t1(prim_gauss,1)


atomlabel='Si'
gbs_path ='/Users/believe/Google\ Drive/unt2/xintegrals/basis/sto3g.gbs'
subblocks = [1]
prim_gauss = fitsubs.get_prim_gauss(atomlabel,gbs_path,subblocks)
#plot_contracted_gauss(prim_gauss,0)
fitsubs.plot_contracted_gauss_t1(prim_gauss,1)
 
fullpath='/Users/believe/Google Drive/unt2/xintegrals/' + 'test.xml'
rmesh, wfc_mat = fitsubs.read_qe_radial_function(fullpath)

#fitting s
#select s gaussian
l = 0
l_index = prim_gauss[:,0] == l
prim_gauss_cart = prim_gauss[l_index,:] #these are not normalized
l_phi   = wfc_mat[0]
print prim_gauss_cart
prim_gauss_sph = fitsubs.prim_spherical_gaussian(prim_gauss_cart)
print prim_gauss_sph
#>>> print prim_gauss_sph
#[[ 0.          1.47874062 -0.74405775]
# [ 0.          0.41256488  0.29340298]
# [ 0.          0.1614751   0.57946728]]
#>>> print prim_gauss_cart
#[[ 0.          1.47874062 -0.21962037]
# [ 0.          0.41256488  0.22559543]
# [ 0.          0.1614751   0.90039843]]


def residuals(l_gauss_coeff,l_gauss_exp,l_phi,r,l):
    err = l_phi/r/r**l 
    counter = 0
    for counter in range(0,l_gauss_exp.size):
     err = err - l_gauss_coeff[counter]*exp(-l_gauss_exp[counter]*r**2)
    return err

from scipy.optimize import leastsq
l_gauss_coeff0 = prim_gauss_cart[:,2]
l_gauss_exp    = prim_gauss_cart[:,1]
r = rmesh
lsq_coeff, other = leastsq(residuals, l_gauss_coeff0, \
                    args=(l_gauss_exp,l_phi,r,l))
print lsq_coeff
#>>> print lsq_coeff
#[ 0.00652299  0.08234451  0.65839353]

l_gauss_coeff0 = prim_gauss_sph[:,2]
l_gauss_exp    = prim_gauss_sph[:,1]
lsq_coeff, other = leastsq(residuals, l_gauss_coeff0, \
                    args=(l_gauss_exp,l_phi,r,l))
print lsq_coeff
#>>> print lsq_coeff
#[ 0.00652299  0.08234451  0.65839353]


#import matlibplot.pyplot as plt
plt.figure()
yyt = np.zeros(rmesh.shape)
igauss = 0
for ii in lsq_coeff:
 yy = ii*np.exp(-l_gauss_exp[igauss]*rmesh**2)
 plt.plot(rmesh,yy,label='ngauss='+str(igauss))
 igauss = igauss + 1
 yyt = yyt +yy
plt.plot(rmesh,yyt,label='ngauss=total')
plt.plot(rmesh,l_phi/rmesh,label='ref:phi(r)/r',ls='dotted')
plt.xlabel('radius (Bohr)')
plt.ylabel(' wfc (Bohr^-3)')
plt.xlim([0,10])
plt.legend()
plt.title('Optimizing only the coefficients')
plt.show()

from scipy.optimize import leastsq
#all exponents, then all coefficients
#l_gauss2_init    = concatenate((l_gauss[:,1],l_gauss[:,2]))
#l_gauss2_init    = concatenate((lsq_coeff,l_gauss[:,2]))
l_gauss2_init    = concatenate((prim_gauss_sph[:,1],prim_gauss_sph[:,2]))
r = rmesh
lsq_coeff2, other = leastsq(residuals_2, l_gauss2_init, \
                    args=(l_phi,r,l))
print lsq_coeff2
#>>> [ 0.21974211  0.21974935  0.08625912 -0.15468558  0.73896893  0.16276638]

#when starting with the u
#>>> print lsq_coeff2
#[ 0.06798315  0.20215311  2.00258805  0.09379241  0.64642417  0.00704821]

plt.figure()
yyt = np.zeros(rmesh.shape)
igauss = 0
lsq_solution = reshape(lsq_coeff2,(2,3)).T
for ii in lsq_solution:
 yy = ii[1]*np.exp(-ii[0]*rmesh**2)
 plt.plot(rmesh,yy,label='ngauss='+str(igauss))
 igauss = igauss + 1
 yyt = yyt +yy

plt.plot(rmesh,yyt,label='ngauss=total')
plt.plot(rmesh,l_phi/rmesh,label='ref:phi(r)/r',ls='dotted')
plt.xlabel('radius (Bohr)')
plt.ylabel(' wfc (Bohr^-3)')
plt.xlim([0,10])
plt.legend()
plt.title('Optimizing both exponents and coefficients')
plt.show()                    

#type of orbital
l = 0
exps_coeffs = lsq_solution

##################################################
#fitting the p to sum of  gaussians
#Mon Jul 22 17:29:10 CDT 2013
atomlabel='Si'
gbs_path ='/Users/believe/Google\ Drive/unt2/xintegrals/basis/sto3g.gbs'
subblocks = [3] #Starts from 1, not zero, i.e. 1 is the first block
prim_gauss = fitsubs.get_prim_gauss(atomlabel,gbs_path,subblocks)

fullpath='/Users/believe/Google Drive/unt2/xintegrals/' + 'test.xml'
rmesh, wfc_mat = fitsubs.read_qe_radial_function(fullpath)

#fitting p (l=1)
l = 1
l_index = prim_gauss[:,0] == l
prim_gauss_cart = prim_gauss[l_index,:] #these are not normalized
l_phi   = wfc_mat[l]    #careful. this is set manually
print prim_gauss_cart
prim_gauss_sph = fitsubs.prim_spherical_gaussian(prim_gauss_cart)
print prim_gauss_sph
#print prim_gauss_cart
#>>> [[ 1.          1.47874062  0.0105876 ]
# [ 1.          0.41256488  0.59516701]
# [ 1.          0.1614751   0.46200101]]
#>>> print prim_gauss_sph
#[[ 1.          1.47874062  0.05036712]
# [ 1.          0.41256488  0.57410134]
# [ 1.          0.1614751   0.13796193]]
opt_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph(prim_gauss_sph,l_phi,rmesh,l)
#print opt_prim_gauss_sph
#[[ 0.20089801  0.23208065]
# [ 0.59253559  0.23196779]
# [ 0.06062894  0.03651938]]

title = 'optimizing both exponents and coefficients'
fitsubs.plot_prim_gauss_sph(opt_prim_gauss_sph,rmesh,l,l_phi,title)

fullpath = '/Users/believe/Google Drive/unt2/xintegrals/'
atomlabel = 'Si'
l=1
fitsubs.print_opt_prim_gauss_cart(opt_prim_gauss_sph,fullpath,atomlabel,l)

#execfile('Si_p.py')

##################################################
#Fitting Zn to Gaussians
##################################################
#Tue Aug 27 12:25:11 CDT 2013

/Users/believe/gdrive/unt2/xintegrals/test.xml
/Users/believe/Dropbox/dailies/2013/07_july/fitsubs.py

import sys
sys.path.append('/Users/believe/Dropbox/dailies/2013/07_july')
import fitsubs

atomlabel='Zn'
#gbs_path ='/Users/believe/gdrive/UNT2/xintegrals/basis/luis_crenbs.gbs'
gbs_path ='/Users/believe/gdrive/UNT2/xintegrals/basis/luis_sto3g.gbs'
#Fitting s
subblocks = [4] #Starts from 1, not zero, i.e. 1 is the first block. This is usually to
                  #pick the valence from the core+valence electrons
prim_gauss = fitsubs.get_prim_gauss(atomlabel,gbs_path,subblocks)

#Treat your UPF file to comply with the xml parser
#1.Works with format UPF 2. For UPF version 1 use the converter
#  $QE/upftools/upf2upf2.x 
#2.Delete any line that has "&" as first non blank character. It confuses they python xml parser. 


fullpath='/Users/believe/gdrive/unt2/xintegrals/O1Zn1_ICSD_26170/qe/Zn.pbe-van_ak_converted.UPF' 
rmesh, wfc_mat = fitsubs.read_qe_radial_function(fullpath)

#fitting s (l=0)
l = 0
l_index = prim_gauss[:,0] == l
prim_gauss_cart = prim_gauss[l_index,:] #these are not normalized
#
#careful. this is set manually. it should correspond to the ordering of the wfc UPF
import matplotlib.pyplot as plt
import os
os.system('clear')
l_phi   = wfc_mat[1]    
title='cartesian gaussians'
fitsubs.plot_prim_gauss(prim_gauss_cart,rmesh,l,l_phi,title)
plt.ylim((-0.5, 1.5))
plt.draw()

prim_gauss_sph = fitsubs.prim_spherical_gaussian(prim_gauss_cart)
title='spherical gaussians'
fitsubs.plot_prim_gauss(prim_gauss_sph,rmesh,l,l_phi,title)
plt.ylim((-0.5, 1.5))
plt.draw()

opt1_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph_coeffs_only(prim_gauss_sph,l_phi,rmesh,l)
title='optimized gaussians: coefficients'
fitsubs.plot_prim_gauss_sph(opt1_prim_gauss_sph,rmesh,l,l_phi,title)
plt.ylim((-0.5, 1.5))
plt.draw()

opt2_prim_gauss_sph = np.concatenate(( np.zeros((3,1)),opt1_prim_gauss_sph ),axis=1)
opt3_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph(opt2_prim_gauss_sph,l_phi,rmesh,l)
title='optimized gaussians: exponents and coefficients from opt coeffs'
fitsubs.plot_prim_gauss_sph(opt3_prim_gauss_sph,rmesh,l,l_phi,title)
plt.ylim((-0.5, 1.5))
plt.draw()

opt_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph(prim_gauss_sph,l_phi,rmesh,l)
title='optimized gaussians: exponents and coefficients'
fitsubs.plot_prim_gauss_sph(opt_prim_gauss_sph,rmesh,l,l_phi,title)
plt.ylim((-0.5, 1.5))
plt.draw()
>>> print prim_gauss_cart
[[ 0.          0.88971389 -0.30884412]
 [ 0.          0.32836038  0.01960641]
 [ 0.          0.14500741  1.13103444]]
>>> print prim_gauss_sph
[[ 0.          0.88971389 -0.7148125 ]
 [ 0.          0.32836038  0.02148702]
 [ 0.          0.14500741  0.67148063]]
>>> print opt1_prim_gauss_sph
[[ 0.88971389 -0.28537333]
 [ 0.32836038 -0.33068241]
 [ 0.14500741  0.77084979]]
>>> print opt_prim_gauss_sph
[[  2.03523495 -28.51986028]
 [  2.05221347  28.17919243]
 [  0.10318553   0.48100614]]
print prim_gauss_cart
print prim_gauss_sph
print opt1_prim_gauss_sph
print opt_prim_gauss_sph

#Retrying to sto-3g
import sys
sys.path.append('/Users/believe/Dropbox/dailies/2013/07_july')
import fitsubs
#atomlabel='Zn'
atomlabel='Cd'
gbs_path ='/Users/believe/gdrive/UNT2/xintegrals/basis/luis_sto3g.gbs'
#Fitting s
subblocks = [2, 5] #Starts from 1, not zero, i.e. 1 is the first block. This is usually to
                  #pick the valence from the core+valence electrons
prim_gauss = fitsubs.get_prim_gauss(atomlabel,gbs_path,subblocks)
#fullpath='/Users/believe/gdrive/unt2/xintegrals/O1Zn1_ICSD_26170/qe/Zn.pbe-van_ak_converted.UPF' 
#rmesh, wfc_mat = fitsubs.read_qe_radial_function(fullpath)
##fitting s (l=0)
l = 0
l_index = prim_gauss[:,0] == l
prim_gauss_cart = prim_gauss[l_index,:] #these are not normalized
#
#careful. this is set manually. it should correspond to the ordering of the wfc UPF
l_phi   = wfc_mat[1]    
print prim_gauss_cart
prim_gauss_sph = fitsubs.prim_spherical_gaussian(prim_gauss_cart)
print prim_gauss_sph
opt_prim_gauss_sph = fitsubs.optimize_prim_gauss_sph(prim_gauss_sph,l_phi,rmesh,l)
print opt_prim_gauss_sph
title = 'optimizing both exponents and coefficients'
fitsubs.plot_prim_gauss_sph(opt_prim_gauss_sph,rmesh,l,l_phi,title)




fullpath = '/Users/believe/Google Drive/unt2/xintegrals/'
atomlabel = 'Si'
l=1
fitsubs.print_opt_prim_gauss_cart(opt_prim_gauss_sph,fullpath,atomlabel,l)


atomlabel='O'
gbs_path ='/Users/believe/gdrive/UNT2/xintegrals/basis/sto3g.gbs'
subblocks = [2] #Starts from 1, not zero, i.e. 1 is the first block
prim_gauss = fitsubs.get_prim_gauss(atomlabel,gbs_path,subblocks)


