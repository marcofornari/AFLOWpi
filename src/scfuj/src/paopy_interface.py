# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOWpi software.
#
#  AFLOWpi is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

import AFLOWpi.prep
import os
import numpy as np 

def _run_paopy(oneCalc,ID,acbn0=False,exec_prefix=""):
    paopy_path = os.path.join(AFLOWpi.__path__[0],'PAOpy/src','test.py')

    if exec_prefix=="":

        if acbn0:
            execPrefix = AFLOWpi.prep._ConfigSectionMap('run','exec_prefix_serial')
        else:
            execPrefix = AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

    else:
        execPrefix=exec_prefix

    paopy_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_PAOpy.out'%ID)

    paopy_input = 'inputfile.py'
    try:
        command = '%s python %s %s > %s' % (execPrefix,paopy_path,paopy_input,paopy_output)
        print command
        os.system(command)
    except Exception,e:
        print e



    


def paopy_header_wrapper(calcs,shift_type=1,shift='auto',thresh=0.90,tb_kp_mult=4):
    content = """AFLOWpi.scfuj._add_paopy_header(oneCalc,ID,shift_type=%s,shift='auto',thresh=%s,tb_kp_mult=%s)""" % \
        (shift_type,thresh,tb_kp_mult)
        
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',content)

def paopy_dos_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_dos(oneCalc,ID)""")

def paopy_pdos_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_pdos(oneCalc,ID)""")

def paopy_bands_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_bands(oneCalc,ID)""")

def paopy_transport_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_transport(oneCalc,ID)""")
    AFLOWpi.prep.addToAll_(calcs,'POSTPROCESSING','AFLOWpi.scfuj._rename_boltz_files(oneCalc,ID)')
def paopy_optical_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_optical(oneCalc,ID)""")
    AFLOWpi.prep.addToAll_(calcs,'POSTPROCESSING','AFLOWpi.scfuj._rename_boltz_files(oneCalc,ID)')

def paopy_acbn0_wrapper(calcs):
    pass


def _add_paopy_header(oneCalc,ID,shift_type=1,shift='auto',thresh=0.90,tb_kp_mult=4,acbn0=False,ovp=False):
    
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.py')
    ibrav=int(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ibrav'])

    nk1,nk2,nk3 = AFLOWpi.scfuj._mult_kgrid(oneCalc,mult=tb_kp_mult)


    with open(paopy_input,'w') as ifo:
        ifo.write('fpath = "%s_TB.save"\n'%ID)
        if shift=="auto":
            ifo.write('shift="auto"\n')
        else:
            ifo.write('shift=%s\n'%shift)
        ifo.write('shift_type = %s\n'%shift_type)
        ifo.write('ibrav = %s\n'%ibrav)
        ifo.write('pthr = %s\n'%thresh)
        ifo.write('verbose = True\n')
        if float(tb_kp_mult)!=1.0:
            ifo.write('double_grid = True\n')
            ifo.write('nfft1 = %s\n'%nk1)
            ifo.write('nfft2 = %s\n'%nk2)
            ifo.write('nfft3 = %s\n'%nk3)
        if acbn0==True:
            ifo.write('write2file = True\n')
            ifo.write('write_binary = True\n')
        if ovp==True:
            ifo.write('non_ortho = True\n')
def _add_paopy_dos(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.py')
    with open(paopy_input,'a') as ifo:
        ifo.write('do_dos = True\n')
        ifo.write('delta = 0.1\n')
        ifo.write('emin = -12.0\n')
        ifo.write('emax =  12.0\n')
def _add_paopy_pdos(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.py')
    with open(paopy_input,'a') as ifo:
        ifo.write('do_pdos = True\n')
        ifo.write('delta = 0.1\n')
        ifo.write('emin = -12.0\n')
        ifo.write('emax =  12.0\n')

def _add_paopy_bands(oneCalc,ID,nk=1000):
    dk = AFLOWpi.retr._getPath_nk2dk(nk, oneCalc,ID=ID)
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.py')
    with open(paopy_input,'a') as ifo:
        ifo.write('do_bands = True\n')
        ifo.write('nk = %s\n'%nk)

def _add_paopy_transport(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.py')
    with open(paopy_input,'a') as ifo:
        ifo.write('Boltzmann = True\n')



def _add_paopy_optical(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.py')
    with open(paopy_input,'a') as ifo:
        ifo.write('epsilon = True\n')
    





def _mult_kgrid(oneCalc,mult=5.0):

    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    kpt_str  = inputDict['K_POINTS']['__content__']    
    k_grid = [int(np.ceil(float(x)*mult)) for x in kpt_str.split()[:3]]

    return k_grid[0],k_grid[1],k_grid[2]
   
def _rename_boltz_files(oneCalc,ID):
    nspin = AFLOWpi.scfuj.chkSpinCalc(oneCalc,ID=ID)
    temperature='300'
   
    conv_dict={}
    if nspin!=1:
        conv_dict['Seebeck_0.dat']  = '%s_PAOpy_seebeck_up_%sK.dat'%(ID,temperature)           
        conv_dict['sigma_0.dat']    = '%s_PAOpy_cond_up_%sK.dat'%(ID,temperature)     
        conv_dict['kappa_0.dat']    = '%s_PAOpy_kappa_up_%sK.dat'%(ID,temperature)             
        conv_dict['epsr_0.dat']     = '%s_PAOpy_epsilon_up_real.dat'%ID                        
        conv_dict['epsi_0.dat']     = '%s_PAOpy_epsilon_up_imag.dat'%ID                        

        conv_dict['Seebeck_1.dat']  = '%s_PAOpy_seebeck_down_%sK.dat'%(ID,temperature)           
        conv_dict['sigma_1.dat']    = '%s_PAOpy_cond_down_%sK.dat'%(ID,temperature)     
        conv_dict['kappa_1.dat']    = '%s_PAOpy_kappa_down_%sK.dat'%(ID,temperature)             
        conv_dict['epsr_1.dat']     = '%s_PAOpy_epsilon_down_real.dat'%ID                        
        conv_dict['epsi_1.dat']     = '%s_PAOpy_epsilon_down_imag.dat'%ID                        
    else:
        conv_dict['Seebeck_0.dat']  = '%s_PAOpy_seebeck_%sK.dat'%(ID,temperature)           
        conv_dict['sigma_0.dat']    = '%s_PAOpy_cond_%sK.dat'%(ID,temperature)     
        conv_dict['kappa_0.dat']    = '%s_PAOpy_kappa_%sK.dat'%(ID,temperature)             
        conv_dict['epsr_0.dat']     = '%s_PAOpy_epsilon_real.dat'%ID                        
        conv_dict['epsi_0.dat']     = '%s_PAOpy_epsilon_imag.dat'%ID                        

    for old,new in conv_dict.iteritems():
        old_split = old.split('_')
        try:
            xx = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],old_split[0]+'_xx_'+old_split[1])
            yy = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],old_split[0]+'_yy_'+old_split[1])
            zz = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],old_split[0]+'_zz_'+old_split[1])
            print xx
            xx_arr = np.loadtxt(xx)
            yy_arr = np.loadtxt(yy)
            zz_arr = np.loadtxt(zz)

            comb_arr = np.concatenate((xx_arr[:,np.newaxis,0],xx_arr[:,np.newaxis,1],
                                       yy_arr[:,np.newaxis,1],zz_arr[:,np.newaxis,1] ),axis=1)
            new_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],new)

            np.savetxt(new_path,comb_arr)

        except Exception,e: 
            try:
                os.rename(old,new)
            except Exception,e: print e
