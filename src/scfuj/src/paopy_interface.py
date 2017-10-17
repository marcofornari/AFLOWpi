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
import re
import itertools

def _run_paopy(oneCalc,ID,acbn0=False,exec_prefix=""):
    paopy_path = os.path.join(AFLOWpi.__path__[0],'PAOpy/src','main.py')

    if exec_prefix=="":

        if acbn0:
            execPrefix = AFLOWpi.prep._ConfigSectionMap('run','exec_prefix_serial')
        else:
            execPrefix = AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

    else:
        execPrefix=exec_prefix

    paopy_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_PAOpy.out'%ID)

    paopy_input = 'inputfile.xml'
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

def paopy_spin_Hall_wrapper(calcs,spin_texture=False):
    command = """AFLOWpi.scfuj._add_paopy_spin_Hall(oneCalc,ID,spin_texture=%s)"""%spin_texture
    for ID,oneCalc in calcs.iteritems():
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.py'%ID),'r') as ifo:
            input_text = ifo.read()
        input_text = re.sub('###_PAOFLOW_SPECIAL_###',command,input_text)
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.py'%ID),'w') as ofo:
            ofo.write(input_text)

def paopy_Berry_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_Berry(oneCalc,ID)""")

def paopy_dos_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_dos(oneCalc,ID)""")

def paopy_pdos_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_pdos(oneCalc,ID)""")

def paopy_bands_wrapper(calcs,band_topology=True,fermi_surface=False):
    AFLOWpi.prep.addToAll_(calcs,
                           'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_bands(oneCalc,ID,topology=%s,fermi_surface=%s)"""%(band_topology,fermi_surface))

def paopy_transport_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_transport(oneCalc,ID)""")
    AFLOWpi.prep.addToAll_(calcs,'POSTPROCESSING','AFLOWpi.scfuj._rename_boltz_files(oneCalc,ID)')
def paopy_optical_wrapper(calcs):
    AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""AFLOWpi.scfuj._add_paopy_optical(oneCalc,ID)""")
    AFLOWpi.prep.addToAll_(calcs,'POSTPROCESSING','AFLOWpi.scfuj._rename_boltz_files(oneCalc,ID)')

def paopy_acbn0_wrapper(calcs):
    pass





def _add_paopy_xml(filename,var_name,var_type,var_val,array=False):  

    with open(filename,'r') as ifo:                                                                                   
        lines=ifo.readlines()
                                                                                         
    var_size=1
    if array:
        var_size = len(var_val)
        var_val = ' '.join(map(str,var_val))
        lines[-1]='\t<%s><a type="%s" size="%s">%s</a></%s>'%(var_name,var_type,var_size,var_val,var_name)            
    else:
        lines[-1]='\t<%s type="%s" size="%s">%s</%s>'%(var_name,var_type,var_size,var_val,var_name)            
    lines.extend(['\n</root>'])                                                                                       
    outstr = ''.join(lines)                                                                                           
    with open(filename,'w') as ofo:                                                                                   
        ofo.write(outstr)


def _add_paopy_header(oneCalc,ID,shift_type=1,shift='auto',thresh=0.90,tb_kp_mult=4,acbn0=False,ovp=False,smearing='gauss'):
    
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    ibrav=int(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ibrav'])

    nk1,nk2,nk3 = AFLOWpi.scfuj._mult_kgrid(oneCalc,mult=tb_kp_mult)

    blank = '''<?xml version="1.0"?>
<root>
</root>'''

    with open(paopy_input,'w') as ifo:
        ifo.write(blank)
    
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'fpath','character','%s_TB.save'%ID)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'smearing','character',smearing)

    if shift=="auto":
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'shift','character',"auto")
    else:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'shift','decimal',shift)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'shift_type','int',shift_type)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ibrav','int',ibrav)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'pthr','decimal',thresh)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'verbose','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'npool','int',1)
    if float(tb_kp_mult)!=1.0:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'double_grid','logical','T')
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nfft1','int',nk1)
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nfft2','int',nk2)
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nfft3','int',nk3)
    if acbn0==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'write2file','logical','T')
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'write_binary','logical','T')
    if ovp==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'non_ortho','logical','T')

def _add_paopy_dos(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'do_dos','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'delta','decimal',0.1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emin','decimal',-12.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emax','decimal',12.0)

def _add_paopy_pdos(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'do_pdos','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'delta','decimal',0.1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emin','decimal',-12.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emax','decimal',12.0)
    
def _add_paopy_bands(oneCalc,ID,nk=1000,topology=True,fermi_surface=False):

    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'do_bands','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nk','int',nk)
    if topology==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'band_topology','logical','T')
    if fermi_surface==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'fermisurf','logical','T')



def _add_paopy_transport(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'Boltzmann','logical','T')

def _add_paopy_optical(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'epsilon','logical','T')

    
def _add_paopy_spin_Hall(oneCalc,ID,spin_texture=False):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')

    nl,sh = AFLOWpi.scfuj._get_spin_ordering(oneCalc,ID)

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'sh','int',sh,array=True)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nl','int',nl,array=True)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'spin_Hall','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'delta','decimal',0.1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'eminSH','decimal',-12.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emaxSH','decimal',12.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ac_cond_spin','logical','T')

    if spin_texture:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'spintexture','logical','T')

def _add_paopy_Berry(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'Berry','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'delta','decimal',0.1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'eminAH','decimal',-12.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emaxAH','decimal',12.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ac_cond_Berry','logical','T')



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




def _get_spin_ordering(oneCalc,ID):
    in_str = AFLOWpi.retr._getOutputString(oneCalc,ID+'_pdos')

    rename_info_re = re.compile(r'state #\s*(\d*): atom\s*(\d+)\s*\(\s*(\S*)\s*\).*wfc\s*(\d+).*l=(\d+).*m_j=([\s-][.\d]+).*\n')

    res = rename_info_re.findall(in_str)

    l_list= []
    for i in res:
         l_list.append(int(i[4]))

    grouped_l = [list(g) for k, g in itertools.groupby(l_list)] 
    sh = []
    nl = []

    for i in xrange(len(grouped_l)):
        sh.append(grouped_l[i][0])
        if grouped_l[i][0]==0:
            nl.append(len(grouped_l[i])/2)
        if grouped_l[i][0]==1:
            nl.append(len(grouped_l[i])/6)
        if grouped_l[i][0]==2:
            nl.append(len(grouped_l[i])/10)
        if grouped_l[i][0]==3:
            nl.append(len(grouped_l[i])/14)


    return nl,sh
