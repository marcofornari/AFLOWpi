#By Luis Agapito at Marco Buongiorno Nardelli's UNT group
#2013
##################################################

import numpy as np

##################################################
def print_opt_prim_gauss_cart(prim_gauss_sph,fullpath,atomlabel,l):
 #creating a dictionary of primitive cartesian gaussians from the fitted spherical gaussians
 import scipy.misc as spm
 # spherical to cartesian gaussians
 #l,m)  : [ (Nlm , [( c, lx,ly,lz)])
 s2c = {
 (0, 0) :  (1/4 , [( 1, 0, 0, 0)] ) ,
 (1,-1) :  (1/4 , [( 1, 0, 1, 0)] ) ,
 (1, 0) :  (1/4 , [( 1, 0, 0, 1)] ) ,
 (1, 1) :  (1/4 , [( 1, 1, 0, 0)] ) ,
 (2,-2) :  (1/4 , [( 1, 1, 1, 0)] ) ,
 (2,-1) :  (1/4 , [( 1, 0, 1, 1)] ) ,
 (2, 0) :  (  3 , [( 2, 0, 0, 2),
                   (-1, 2, 0, 0),
                   (-1, 0, 2, 0)] ) ,
 (2, 1) :  (1/4 , [( 1, 1, 0, 1)] ) ,
 (2, 2) :  (  1 , [( 1, 2, 0, 0),
                   (-1, 0, 2, 0)] ) }
 
 PCG = {}
 for m in range(-l,l+1):
  Nlm, mytuple = s2c[(l,m)]
  g_tuples = []
  for counter,(d,lx,ly,lz) in enumerate(mytuple):
   factor_lm = 1/4*np.sqrt(spm.factorial2(2*l+1)/np.pi/Nlm)
   for prim_g in prim_gauss_sph: 
    print(prim_g)
    g_tuples = g_tuples + [(lx, ly, lz, d*prim_g[1]*factor_lm,prim_g[0])]
   dict_field= {(l,m) : g_tuples} 
   PCG.update(dict_field)

 llabel = { 0 : 's', 1 : 'p', 2 : 'd' }
 fullname = fullpath+'/'+atomlabel+'.'+llabel[l]+'.py'
 f = open(fullname,'w+')
 #f.write(atomlabel+'.'+llabel[l]+ '=' + repr(PCG) + '\n' )
 f.write('cgto=' + repr(PCG) + '\n' )
 f.close()

##################################################
def plot_prim_gauss(prim_gauss,rmesh,l,l_phi,title):
 import matplotlib.pyplot as plt
 plt.figure(figsize=(5,4))
 yyt = np.zeros(rmesh.shape)
 igauss = 0
 for ii in prim_gauss:
  yy = ii[2]*np.exp(-ii[1]*rmesh**2)*rmesh**(l)
  plt.plot(rmesh,yy,label='L='+str(l)+' ngauss='+str(igauss),)
  igauss = igauss + 1
  yyt = yyt +yy
 
 plt.plot(rmesh,yyt,label='ngauss=total')
 plt.plot(rmesh,l_phi/rmesh,label='ref:phi(r)/r',ls='dotted',color='k')
 plt.xlabel('radius (Bohr)')
 plt.ylabel(' wfc (Bohr^-3)')
 plt.xlim([0,10])
 plt.legend()
 plt.title(title)
 plt.draw()                    
 plt.show(block=False)                    

##################################################
def plot_prim_gauss_sph(prim_gauss_sph,rmesh,l,l_phi,title):
 import matplotlib.pyplot as plt
 plt.figure(figsize=(5,4))
 yyt = np.zeros(rmesh.shape)
 igauss = 0
 for ii in prim_gauss_sph:
  yy = ii[1]*np.exp(-ii[0]*rmesh**2)*rmesh**(l)
  plt.plot(rmesh,yy,label='L='+str(l)+' ngauss='+str(igauss),)
  igauss = igauss + 1
  yyt = yyt +yy
 
 plt.plot(rmesh,yyt,label='ngauss=total')
 plt.plot(rmesh,l_phi/rmesh,label='ref:phi(r)/r',ls='dotted',color='k')
 plt.xlabel('radius (Bohr)')
 plt.ylabel(' wfc (Bohr^-3)')
 plt.xlim([0,10])
 plt.legend()
 plt.title(title)
 plt.draw()
 plt.show(block=False)                    

##################################################
def optimize_prim_gauss_sph_coeffs_only(prim_gauss_sph,l_phi,rmesh,l):
 from scipy.optimize import leastsq
 l_coeffs_init     = prim_gauss_sph[:,2]
 l_exps            = prim_gauss_sph[:,1]
 print(("initial coefficients:",l_coeffs_init))
 print(("initial exps        :",l_exps))
 print(("initial residuals   :",sum(residuals_1(l_coeffs_init,l_exps,l_phi,rmesh,l)**2)))
 
#lsq_coeffs, other = leastsq(residuals_1, l_coeffs_init, args=(l_exps,l_phi,rmesh,l))
 lsq_coeffs, cov_x,infodict,mesg,ier = leastsq(residuals_1, l_coeffs_init, args=(l_exps,l_phi,rmesh,l),full_output=1,ftol=1e-9,xtol=1e-9)

 aux = np.concatenate((l_exps,lsq_coeffs))
 print(("final coefficients  :",lsq_coeffs))
 print(("final exps          :",l_exps))
 print(("final residuals     :",sum(residuals_1(lsq_coeffs,l_exps,l_phi,rmesh,l)**2)))
 return np.reshape(aux,(2,-1)).T

##################################################
def optimize_prim_gauss_sph(prim_gauss_sph,l_phi,rmesh,l):
 from scipy.optimize import leastsq
 print(("initial coefficients:",prim_gauss_sph[:,2]))
 print(("initial exps        :",prim_gauss_sph[:,1]))
 l_gauss2_init    = np.concatenate((prim_gauss_sph[:,1],prim_gauss_sph[:,2]))
 print(("initial residuals   :",sum(residuals_2(l_gauss2_init,l_phi,rmesh,l)**2)))

 #lsq_coeff2, other = leastsq(residuals_2, l_gauss2_init, args=(l_phi,rmesh,l))
 lsq_coeff2, cov_x,infodict,mesg,ier = leastsq(residuals_2, l_gauss2_init, args=(l_phi,rmesh,l),full_output=1,ftol=1e-9,xtol=1e-9)
 print(("number of calls:",infodict['nfev']))
 if ier>4:
    print(("ier=",ier,"  Error message:",mesg))
 solut = np.reshape(lsq_coeff2,(2,-1)).T


 print(("final coefficients  :",solut[:,1]))
 print(("final exps          :",solut[:,0]))
 print(("final residuals     :",sum(residuals_2(lsq_coeff2,l_phi,rmesh,l)**2)))
 return solut 
##################################################
#residuals_1 optimizes coefficients but not exponents
def residuals_1(l_coeffs,l_exps_fixed,l_phi,r,l):
    #err = l_phi/r/r**l 
    err = l_phi 
    counter = 0
    nprimgauss = int(l_coeffs.size);
    for counter in range(0,nprimgauss):
     err = err - r**(l+1)*l_coeffs[counter]*np.exp(-l_exps_fixed[counter]*r**2)
    return err

##################################################
#residuals_2 also optimizes the exponents
def residuals_2(l_gauss2,l_phi,r,l):
    #err = l_phi/r/r**l 
    err = l_phi 
    counter = 0
    nprimgauss = int(l_gauss2.size/2);
    for counter in range(0,nprimgauss):
     err = err - r**(l+1)*l_gauss2[counter+nprimgauss]*np.exp(-l_gauss2[counter]*r**2)
    return err

##################################################
def prim_spherical_gaussian(prim_gauss_cart):
    import numpy as np
    import copy
    prim_gauss_spherical = copy.copy(prim_gauss_cart) # l, zeta, c
    
#for d it takes the normalization factor of the dxy , which is = dyz = dxz != dx2-y2 ,etc
    #factors = i!j!k!/(2i)!/(2j)!/(2k)!
    factorials = { 0 : 1, 1 : 1/2, 2: 1/4 }
    print((factorials[2]))
    c2s_factor = { 0 : 3.54490770181103, 1 : 2.04665341589298, 2 : 0.915291232863769 }
    #>>> c2s[(0,0,0)]()
    #3.54490770181103*s*exp(-r**2)
    #>>> c2s[(1,0,0)]()
    #2.04665341589298*px*r*exp(-r**2)
    #>>> c2s[(1,1,0)]()
    #0.915291232863769*dxy*r**2*exp(-r**2) 
    for counter,(l,alphai,ci) in enumerate(prim_gauss_cart):
     prim_gauss_spherical[counter,2] = (2*alphai/np.pi)**(3/4)*np.sqrt( (8*alphai)**(l)*factorials[l] )\
                                       * c2s_factor[l] * ci
    return prim_gauss_spherical  #these contain the normalization factor of the cartesian gaussians


##################################################
def read_qe_radial_function(fullpath):
 import numpy as np
 import re
 import xml.etree.ElementTree as et
 import sys
 #reads xml file
 #http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format
 Bohr2Angs=0.529177249000000
 #fpath='/Users/believe/Google Drive/unt2/xintegrals/'
 #tree = et.parse(fpath+'Si.pbe-n-rrkjus_psl.0.1.UPF')
 #tree = et.parse(fpath+'test.xml')
 tree = et.parse(fullpath)
 root = tree.getroot()
  
 for xx in root.iter('PP_R'):
   nmesh1 = int(xx.attrib['size'])
   xxaux = re.split('\n| ',xx.text)
 rmesh  =np.array(list(map(float,[_f for _f in xxaux if _f]))) #In Bohrs
 nmesh  =len(rmesh) 
 if nmesh1 != nmesh:
   sys.exit('Error deciding the size of the r mesh')

 print(('Number of radial points: %i' % nmesh))
 
 xvar = list(range(1,nmesh+1))
 #plt.plot(xvar,rmesh*Bohr2Angs,'ro')
 #plt.yscale('log')
 #plt.xlabel('mesh counter')
 #plt.ylabel('radius (Angs)')
 #plt.show()
 
 wfc_mat = np.zeros((0,nmesh))
 pswfc = root.find('PP_PSWFC')
 labels = []
  
 for node in pswfc:
   print((node.tag, node.attrib['l'],node.attrib['label']))
   labels.append('l='+node.attrib['l']+'   '+node.attrib['label'])
   xxaux = re.split('\n| ',node.text)
   wfc_aux  =np.array([list(map(float,[_f for _f in xxaux if _f]))])
   wfc_mat = np.concatenate((wfc_mat,wfc_aux))
 
 print(('Number of radial wavefunctions found: %i' % wfc_mat.shape[0]))
 
 return rmesh, wfc_mat #in Bohrs
 
 #for i in range(0,wfc_mat.shape[0]):
 # plt.plot(rmesh*Bohr2Angs,wfc_mat[i]*rmesh,label=labels[i])
 #plt.legend()
 #plt.ylabel('phi(r)*r')  #units are ???
 #plt.xlabel('radius (Angs)')
 #plt.show()

##################################################
def plot_contracted_gauss(prim_gauss,nfig):
 import matplotlib.pyplot as plt
 import numpy as np
 plt.figure(nfig)
 for ns in np.unique(prim_gauss[:,[0]]):
  i,j = np.where(prim_gauss[:,[int(ns)]] == ns)
  wfc = prim_gauss[i,:]
  
  yy = np.zeros(rmesh.shape)
  for ii in wfc:
   yy = yy + (2*ii[1]/np.pi)**(3/4)*ii[2]*np.exp(-ii[1]*rmesh**2)
  plt.plot(rmesh*Bohr2Angs,yy,label='L='+str(int(ns)))
 plt.xlabel('radius (Bohr)')
 plt.ylabel(' wfc (Bohr^3)')
 plt.legend()
 plt.xlim([0,5])
 plt.ylim([0,0.6])
 plt.show(nfig)

##################################################
def plot_contracted_gauss_t1(prim_gauss,nfig):
 import matplotlib.pyplot as plt
 import numpy as np
 plt.figure(nfig)
 sf = 1/1.24**2
 for ns in np.unique(prim_gauss[:,[0]]):
  i,j = np.where(prim_gauss[:,[int(ns)]] == ns)
  wfc = prim_gauss[i,:]
  
  yy = np.zeros(rmesh.shape)
  #it turns out that the gauss table was fitted to \zeta=1.24 already.
  #to match the examples one needs to use \zeta=1. This implies a 
  #division by 1.24^2 (see Levine, p.490)
  for ii in wfc:
   yy = yy + (2*ii[1]*sf/np.pi)**(3/4)*ii[2]*np.exp(-ii[1]*sf*rmesh**2)
  plt.plot(rmesh,yy,label='L='+str(int(ns)))
 plt.xlabel('radius (Bohr)')
 plt.ylabel(' wfc (Bohr^3)')
 plt.legend()
 plt.xlim([0,4])
 plt.ylim([0,0.6])
 plt.show(nfig)

##################################################
def get_prim_gauss(atomlabel,gbs_path,subblocks):
 import subprocess
 import numpy as np
 import re
 import sys
 #subblocks = [3,4]
 #sed_str = 'sed -n "/^-%s/,/^\*\*\*\*/p" %s' % ('Si','/Users/believe/Google\ Drive/unt2/xintegrals/basis/631.gbs')
 sed_str = 'sed -n "/^-%s$/,/^\*\*\*\*$/p" %s' % (atomlabel,gbs_path)
 unix_out = subprocess.check_output(sed_str,shell=True)
 block  = re.split('\n', unix_out)
 print(unix_out)
 print('************')
 prevline = 0
 counter = 1
 bb5=np.zeros((0,3))
 while (prevline+1 < len(block)-2 ): #2 accounts for the last spurious lines in block i.e. ****,\n
  #print prevline+1
  print((' '.join(block[prevline+1].split()).split()))
  shell,nprim,unknown=' '.join(block[prevline+1].split()).split()
  print((shell,nprim,unknown))
  subb= block[prevline+2:prevline+2+int(nprim)]
  bb1=' '.join(subb).replace('D','E')
  bb2=np.array(list(map(float,bb1.split())))
  
  #subb= map(split,block[prevline+2:prevline+2+int(nprim)])
  #print subb
  if (shell == 'S') or (shell == 'P') or (shell == 'D'):
    nshell = 2
  elif shell == 'SP':
    nshell = 3
  elif shell == 'SPD':
    nshell = 4
  else:
   print('case not contemplated')  
  bb3=np.reshape(bb2,(int(nprim),nshell)) 
  #print bb3 
  if any( x in subblocks for x in [counter]):
   #print 'I belong to the chosen subblock',bb3
   #bb4 = np.zeros(shape=(int(nprim)*(nshell-1),2)
   for ii in range(1,nshell):
    bbb1 = bb3[:,[0,ii]]
    if nshell==2:
      if shell=='S':
        ii_shell=1
      elif shell=='P':
        ii_shell=2
      elif shell=='D':
        ii_shell=3
      else:
        print(shell)
        sys.exit('Case not contemplated 2')
    else:
      ii_shell=ii
    bbb2 = np.hstack([np.ones((int(nprim),1))*int(ii_shell-1),bbb1])
    bb5=np.vstack([bb5,bbb2])
 
  prevline = prevline+1 + int(nprim)
  counter = counter +1
 #print 'Printing last bb5',bb5
 return bb5
