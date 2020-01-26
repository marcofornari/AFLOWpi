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
import shutil

def _run_paopy(oneCalc,ID,acbn0=False,exec_prefix=""):
    paopy_path = os.path.join(AFLOWpi.__path__[0],'PAOFLOW/examples/','main.py')

    shutil.copy(paopy_path,oneCalc['_AFLOWPI_FOLDER_'])

    paopy_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'main.py')

    if exec_prefix=="":

        if acbn0:
            execPrefix = AFLOWpi.prep._ConfigSectionMap('run','exec_prefix_serial')
        else:
            execPrefix = AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

    else:
        execPrefix=exec_prefix

    paopy_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_PAOFLOW.out'%ID)

    paopy_input = '%s inputfile.xml'%oneCalc['_AFLOWPI_FOLDER_']
    try:
        command = '%s python -u %s %s > %s' % (execPrefix,paopy_path,paopy_input,paopy_output)
        print(command)
        os.system(command)
    except Exception as e:
        print(e)


def paopy_header_wrapper(calcs,shift_type=1,shift='auto',thresh=0.90,tb_kp_mult=4,smearing=None,emin=-5.0,emax=5.0,ne=1000,symmetrize=False,sym_thr=1.e-6,sym_max_iter=20):
    command="""     AFLOWpi.scfuj._add_paopy_header(oneCalc,ID,shift_type=%s,shift='auto',thresh=%s,tb_kp_mult=%s,smearing=%s,emin=%s,emax=%s,ne=%s,symmetrize=%s,sym_thr=%s,sym_max_iter=%s)""" % (shift_type,thresh,tb_kp_mult,repr(smearing),emin,emax,ne,symmetrize,sym_thr,sym_max_iter)
        
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)

def paopy_spin_Hall_wrapper(calcs,s_tensor,spin_texture=False):
    command="""     AFLOWpi.scfuj._add_paopy_spin_Hall(oneCalc,ID,%s,spin_texture=%s)"""%(repr(s_tensor),spin_texture)
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)


def paopy_Berry_wrapper(calcs,a_tensor):
    command="""     AFLOWpi.scfuj._add_paopy_Berry(oneCalc,ID,%s)"""%repr(a_tensor)
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)

def paopy_dos_wrapper(calcs):
    command="""     AFLOWpi.scfuj._add_paopy_dos(oneCalc,ID)"""
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)

def paopy_pdos_wrapper(calcs):
    command="""     AFLOWpi.scfuj._add_paopy_pdos(oneCalc,ID)"""
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)

def paopy_bands_wrapper(calcs,band_topology=True,fermi_surface=False,ipol=0,jpol=1,spol=2,nk=1000):
    command="""     AFLOWpi.scfuj._add_paopy_bands(oneCalc,ID,topology=%s,fermi_surface=%s,ipol=%s,jpol=%s,spol=%s,nk=%s)"""%(band_topology,fermi_surface,ipol,jpol,spol,nk)
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)

def paopy_transport_wrapper(calcs,t_tensor,t_min,t_max,t_step):
    command="""     AFLOWpi.scfuj._add_paopy_transport(oneCalc,ID,%s,t_min=%s,t_max=%s,t_step=%s)"""%(repr(t_tensor),t_min,t_max,t_step)
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)
    AFLOWpi.prep.addToAll_(calcs,'POSTPROCESSING','AFLOWpi.scfuj._rename_boltz_files(oneCalc,ID)')

def paopy_optical_wrapper(calcs,d_tensor):
    command="""     AFLOWpi.scfuj._add_paopy_optical(oneCalc,ID,%s)"""%repr(d_tensor)
    AFLOWpi.prep.addToAll_(calcs,'PAOFLOW',command)
    AFLOWpi.prep.addToAll_(calcs,'POSTPROCESSING','AFLOWpi.scfuj._rename_boltz_files(oneCalc,ID)')

def paopy_acbn0_wrapper(calcs):
    pass





def _add_paopy_xml(filename,var_name,var_type,var_val,degree=0):  

    with open(filename,'r') as ifo:                                                                 
        lines=ifo.readlines()
                                                                                         
    var_size=1
    if degree==2:
        tensor = '\t<%s>\n'%var_name
        for i in range(len(var_val)):
            tensor+='\t\t<a type="%s" size="%s">%s</a>\n'%(var_type,len(var_val[i]),' '.join(map(str,var_val[i])))
        tensor+="\t</%s>"%var_name
        lines[-1]=tensor        
    elif degree==1:
        var_size = len(var_val)
        var_val = ' '.join(map(str,var_val))
        lines[-1]='\t<%s><a type="%s" size="%s">%s</a></%s>'%(var_name,var_type,var_size,var_val,var_name)   
    else:
        lines[-1]='\t<%s type="%s" size="%s">%s</%s>'%(var_name,var_type,var_size,var_val,var_name)          
    lines.extend(['\n</root>'])                                                                     
    outstr = ''.join(lines)                                                                            
    with open(filename,'w') as ofo:                                                                   
        ofo.write(outstr)


def _add_paopy_header(oneCalc,ID,shift_type=1,shift='auto',thresh=0.90,tb_kp_mult=4,acbn0=False,ovp=False,smearing=None,emin=-5.0,emax=5.0,ne=1000,symmetrize=False,sym_thr=1.e-6,sym_max_iter=20):
    
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    ibrav=oneCalc['_AFLOWPI_ORIG_IBRAV_']
    
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
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'delta','decimal',0.1)

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ne','int',ne)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emin','decimal',emin)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emax','decimal',emax)



    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'double_grid','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nfft1','int',nk1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nfft2','int',nk2)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nfft3','int',nk3)

    if symmetrize:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'symmetrize','logical','T')
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'symm_max_iter','',sym_max_iter)
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'symm_thresh','',sym_thr)


    if acbn0==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'expand_wedge','logical','F')
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'write2file','logical','T')
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'write_binary','logical','T')
    if ovp==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'non_ortho','logical','T')

def _add_paopy_dos(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'do_dos','logical','T')


def _add_paopy_pdos(oneCalc,ID):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'do_pdos','logical','T')

    
def _add_paopy_bands(oneCalc,ID,nk=1000,topology=True,fermi_surface=False,ipol=0,jpol=1,spol=2):

    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'do_bands','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nk','int',nk)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ipol','int',ipol)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'jpol','int',jpol)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'spol','int',spol)
    if topology==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'band_topology','logical','T')
    if fermi_surface==True:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'fermisurf','logical','T')

    if oneCalc['_AFLOWPI_ORIG_IBRAV_']==0:
        HSP,band_path = AFLOWpi.retr._getHighSymPoints(oneCalc,ID)
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'band_path','character',band_path)
        temp_HSP_list=[]
        for k,v in list(HSP.items()):
            tmp = [k]
            tmp.extend(list(map(str,v)))
            temp_HSP_list.append(tmp)

        HSP_ARRAY= np.asarray(temp_HSP_list)

        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'high_sym_points','string',HSP_ARRAY,degree=2)


def _add_paopy_transport(oneCalc,ID,t_tensor,t_min,t_max,t_step):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'Boltzmann','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'tmin','decimal',t_min)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'tmax','decimal',t_max)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'tstep','decimal',t_step)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'t_tensor','int',t_tensor,degree=2)

def _add_paopy_optical(oneCalc,ID,d_tensor):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'epsilon','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'d_tensor','int',d_tensor,degree=2)
    
def _add_paopy_spin_Hall(oneCalc,ID,s_tensor,spin_texture=False):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')

    nl,sh = AFLOWpi.scfuj._get_spin_ordering(oneCalc,ID)

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'sh','int',sh,degree=1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'nl','int',nl,degree=1)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'s_tensor','int',s_tensor,degree=2)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'spin_Hall','logical','T')

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'eminSH','decimal',-5.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emaxSH','decimal',5.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ac_cond_spin','logical','T')

    if spin_texture:
        AFLOWpi.scfuj._add_paopy_xml(paopy_input,'spintexture','logical','T')

def _add_paopy_Berry(oneCalc,ID,a_tensor):
    paopy_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'inputfile.xml')

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'Berry','logical','T')

    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'eminAH','decimal',-5.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'emaxAH','decimal',5.0)
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'ac_cond_Berry','logical','T')
    AFLOWpi.scfuj._add_paopy_xml(paopy_input,'a_tensor','int',a_tensor,degree=2)


def _mult_kgrid(oneCalc,mult=5.0):

    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    kpt_str  = inputDict['K_POINTS']['__content__']    
    try:
        mult[0]
        tmp_kps = kpt_str.split()[:3]
        k_grid = [int(np.ceil(float(tmp_kps[x])*mult[x])) for x in range(len(tmp_kps))]
    except:
        k_grid = [int(np.ceil(float(x)*mult)) for x in kpt_str.split()[:3]]

    return k_grid[0],k_grid[1],k_grid[2]
   
def _rename_boltz_files(oneCalc,ID):
    nspin = AFLOWpi.scfuj.chkSpinCalc(oneCalc,ID=ID)

   
    try:
        test_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'Seebeck_0.dat')
        test_dat  = np.loadtxt(test_file)
        temp = np.unique(test_dat[:,0])
    except Exception as e: 
        return
 
    for T in range(temp.shape[0]):
        conv_dict={}
        if nspin!=1:
            conv_dict['Seebeck_0.dat']  = '%s_PAOFLOW_seebeck_up_%sK.dat'%(ID,int(temp[T]))           
            conv_dict['sigma_0.dat']    = '%s_PAOFLOW_cond_up_%sK.dat'%(ID,int(temp[T]))     
            conv_dict['kappa_0.dat']    = '%s_PAOFLOW_kappa_up_%sK.dat'%(ID,int(temp[T]))             
            conv_dict['epsr_0.dat']     = '%s_PAOFLOW_epsilon_up_real.dat'%ID                        
            conv_dict['epsi_0.dat']     = '%s_PAOFLOW_epsilon_up_imag.dat'%ID                        

            conv_dict['Seebeck_1.dat']  = '%s_PAOFLOW_seebeck_down_%sK.dat'%(ID,int(temp[T]))           
            conv_dict['sigma_1.dat']    = '%s_PAOFLOW_cond_down_%sK.dat'%(ID,int(temp[T]))     
            conv_dict['kappa_1.dat']    = '%s_PAOFLOW_kappa_down_%sK.dat'%(ID,int(temp[T]))             
            conv_dict['epsr_1.dat']     = '%s_PAOFLOW_epsilon_down_real.dat'%ID                        
            conv_dict['epsi_1.dat']     = '%s_PAOFLOW_epsilon_down_imag.dat'%ID                        
        else:
            conv_dict['Seebeck_0.dat']  = '%s_PAOFLOW_seebeck_%sK.dat'%(ID,int(temp[T]))           
            conv_dict['sigma_0.dat']    = '%s_PAOFLOW_cond_%sK.dat'%(ID,int(temp[T]))     
            conv_dict['kappa_0.dat']    = '%s_PAOFLOW_kappa_%sK.dat'%(ID,int(temp[T]))             
            conv_dict['epsr_0.dat']     = '%s_PAOFLOW_epsilon_real.dat'%ID                        
            conv_dict['epsi_0.dat']     = '%s_PAOFLOW_epsilon_imag.dat'%ID                        

        for old,new in list(conv_dict.items()):
            old_split = old.split('_')
            try:
                xx = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],old_split[0]+'_xx_'+old_split[1])
                yy = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],old_split[0]+'_yy_'+old_split[1])
                zz = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],old_split[0]+'_zz_'+old_split[1])
          
                xx_arr = np.loadtxt(xx)
                yy_arr = np.loadtxt(yy)
                zz_arr = np.loadtxt(zz)

                comb_arr = np.concatenate((xx_arr[:,np.newaxis,0],xx_arr[:,np.newaxis,1],
                                           yy_arr[:,np.newaxis,1],zz_arr[:,np.newaxis,1] ),axis=1)
                new_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],new)

                np.savetxt(new_path,comb_arr)

            except Exception as e: 
                try:
                    dat = np.loadtxt(old)
                    temp_dat = dat[np.where(dat[:,0]==temp[T])]
                    np.savetxt(new,temp_dat[:,1:])
                except:
                    pass




def _get_spin_ordering(oneCalc,ID):
    in_str = AFLOWpi.retr._getOutputString(oneCalc,ID+'_pdos')

    rename_info_re = re.compile(r'state #\s*(\d*): atom\s*(\d+)\s*\(\s*(\S*)\s*\).*wfc\s*(\d+).*l=(\d+).*m_j=([\s-][.\d]+).*\n')

    res = rename_info_re.findall(in_str)

    l_list= []
    for i in res:
         l_list.append([int(i[1]),int(i[4])])

    grouped_l = [list(g) for k, g in itertools.groupby(l_list)] 
    sh = []
    nl = []
    print(grouped_l)
    for i in range(len(grouped_l)):
        sh.append(grouped_l[i][0][1])
        if grouped_l[i][0][1]==0:
            nl.append(len(grouped_l[i])/2)
        if grouped_l[i][0][1]==1:
            nl.append(len(grouped_l[i])/6)
        if grouped_l[i][0][1]==2:
            nl.append(len(grouped_l[i])/10)
        if grouped_l[i][0][1]==3:
            nl.append(len(grouped_l[i])/14)

    return nl,sh
