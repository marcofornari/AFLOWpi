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

import AFLOWpi
import os
import re
import __main__
import logging
import shutil
import glob
import numpy
import time
from collections import OrderedDict

class tight_binding:
    def __init__(self,calcs,cond_bands=True,proj_thr=0.95,kp_factor=2.0,proj_sh=5.5,tb_kp_mult=4,exec_prefix="",band_mult=1.0,smearing='gauss',emin=-5.0,emax=5.0,ne=1000,symmetrize=False,sym_thr=1.e-6,sym_max_iter=20):
        self.calcs=calcs
        self.plot=AFLOWpi.prep.tb_plotter(self.calcs)
        self.cond_bands=cond_bands
        self.do_ham=False
        self.step_counter=0
        self.thresh=proj_thr
        self.shift=proj_sh
        self.cond_bands_proj=True
        tetra_nscf=False

        tb_plotter=AFLOWpi.prep.tb_plotter(calcs)


        AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""oneCalc,ID=AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"scf"')""")


        pwx_dir=AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
            symlink=True
        else:
            symlink=False
        pwx_exec_loc = os.path.join(pwx_dir,'pw.x')
        if not os.path.exists(pwx_exec_loc):
            logging.error('ERROR: engine executables not found in %s please check your config file. EXITING' % pwx_dir)
            print(('ERROR: engine executables not found in %s please check your config file EXITING' % pwx_dir))
            raise SystemExit
        AFLOWpi.prep.totree(pwx_exec_loc, calcs,rename=None,symlink=symlink)




        command='''if oneCalc["__execCounter__"]<=%s:
     oneCalc,ID=AFLOWpi.prep._run_tb_ham_prep(__submitNodeName__,oneCalc,ID,kp_factor=%s,band_factor=%s,tetra_nscf=%s)
     oneCalc['__execCounter__']+=1
     AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''%(self.step_counter,kp_factor,band_mult,tetra_nscf)

        AFLOWpi.prep.addToAll_(self.calcs,'RUN',command) 
        self.step_counter+=1

        command='''if oneCalc["__execCounter__"]<=%s:
#PAOFLOW_BLOCK

#END_PAOFLOW_BLOCK

     oneCalc['__execCounter__']+=1
     AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''%(self.step_counter)




        AFLOWpi.prep.addToAll_(self.calcs,'RUN',command)
        self.step_counter+=1

        AFLOWpi.scfuj.paopy_header_wrapper(self.calcs,shift_type=1,shift='auto',
                                           thresh=proj_thr,tb_kp_mult=tb_kp_mult,
                                           smearing=smearing,emin=emin,emax=emax,ne=ne,
                                           symmetrize=symmetrize,sym_thr=sym_thr,
                                           sym_max_iter=sym_max_iter)

        command='''if oneCalc["__execCounter__"]<=%s:
     AFLOWpi.scfuj._run_paopy(oneCalc,ID,exec_prefix="%s")
     AFLOWpi.scfuj.PAOFLOW_DATA_CONV(oneCalc,ID)
     oneCalc['__execCounter__']+=1
     AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''%(self.step_counter,exec_prefix)

        AFLOWpi.prep.addToAll_(self.calcs,'RUN',command)
        self.step_counter+=1


    def shc(self,s_tensor=None,en_range=[0.05,5.05],de=0.05,spin_texture=False):
        '''
        do spin hall conductivity calculation with PAOFLOW

        Keyword Arguments:
              a_tensor (numpy.array): an Nx2 array with the direction to calculate the SHC. 
                                      (default is [[0,1,2],] which is x,y,z direction)
              en_range (array): floats of min and max for energy range
              de (float): energy step
              spin_texture (bool): whether or not to calculate spin texture

        Returns:
              None

        '''             



        if s_tensor is None:
            s_tensor = [[0,0,0],[0,1,0],[0,2,0],[1,0,0],[1,1,0],[1,2,0],[2,0,0],[2,1,0],[2,2,0], \
                            [0,0,1],[0,1,1],[0,2,1],[1,0,1],[1,1,1],[1,2,1],[2,0,1],[2,1,1],[2,2,1], \
                            [0,0,2],[0,1,2],[0,2,2],[1,0,2],[1,1,2],[1,2,2],[2,0,2],[2,1,2],[2,2,2]]
            s_tensor = [[0,1,2],]


        ne=float(en_range[1]-en_range[0])/de

        if self.step_counter==1:
            self.do_ham=True
        else:
            self.do_ham=False

        AFLOWpi.scfuj.paopy_spin_Hall_wrapper(self.calcs,s_tensor,spin_texture=spin_texture)


        calc_type='Spin Hall Conductivity'
        print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                                            AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))

        if spin_texture:
            calc_type='Spin Texture'
            print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))

    def ahc(self,a_tensor=None,en_range=[0.05,5.05],de=0.05):
        '''
        do anomalous hall conductivity calculation with PAOFLOW

        Keyword Arguments:
              a_tensor (numpy.array): an Nx2 array with the direction to calculate the AHC. 
                                      (default is [[0,1],] which is x,y direction)
              en_range (array): floats of min and max for energy range
              de (float): energy step

        Returns:
              None

        '''             


        # Berry curvature
        if a_tensor is None:
            a_tensor = [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]]
            a_tensor = [[0,1],]

        ne=float(en_range[1]-en_range[0])/de

        if self.step_counter==1:
            self.do_ham=True
        else:
            self.do_ham=False
        AFLOWpi.scfuj.paopy_Berry_wrapper(self.calcs,a_tensor)



        calc_type='Anomalous Hall Conductivity'
        print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                                            AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))



    def optical(self,d_tensor=None,en_range=[0.05,5.05],de=0.05):
        '''
        do epsilon calculation with PAOFLOW

        Keyword Arguments:
              a_tensor (numpy.array): an Nx2 array with the direction to calculate epsilon. 
                                      (default is [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]])
              en_range (array): floats of min and max for energy range
              de (float): energy step

        Returns:
              None
        '''             


        if d_tensor is None:
            d_tensor = [[0,0],[1,1],[2,2]]
#            d_tensor = [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]]

        ne=float(en_range[1]-en_range[0])/de

        if self.step_counter==1:
            self.do_ham=True
        else:
            self.do_ham=False
        AFLOWpi.scfuj.paopy_optical_wrapper(self.calcs,d_tensor)

        calc_type='Optical Properties'
        print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                                            AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))





    def transport(self,t_min=300,t_max=300,t_step=1,en_range=[-5.05,5.05],de=0.05,carr_conc=False):
        '''
        Calculate Boltzmann transport properties with PAOFLOW

        Keyword Arguments:
              t_min (float): min of temp range to calculate
              t_max (float): max of temp range to calculate
              t_step (float): temperature step
              en_range (array): floats of min and max for energy range
              de (float): energy step

        Returns:
              None

        '''             

        t_tensor=None

        # Boltzmann transport
        if t_tensor is None:
            t_tensor = [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]]


        ne=float(en_range[1]-en_range[0])/de
        AFLOWpi.scfuj.paopy_transport_wrapper(self.calcs,t_tensor,t_min,t_max,
                                              t_step,carr_conc=carr_conc)

        calc_type='Transport Properties'
        calc_type_color = AFLOWpi.run._colorize_message(calc_type,
                                                        level='DEBUG',
                                                        show_level=False)


        add_step_str = AFLOWpi.run._colorize_message('ADDING TB STEP:  ',
                                                     level='GREEN',
                                                     show_level=False)
        print((add_step_str + calc_type_color))

        return None




    def dos(self,dos_range=[-5.5,5.5],projected=True,de=0.05,fermi_surface=False):
        '''
        do anomalous hall conductivity calculation with PAOFLOW

        Keyword Arguments:
              dos_range (list): energy range around fermi energy to calculate the DOS
              projected (bool): include the projected dos or not
              de (float): energy sampling width
              fermi_surface (bool): whether or not to calculate the fermi surface

        Returns:
              None
        '''       

        k_grid=None
        cond_bands=True

        AFLOWpi.scfuj.paopy_dos_wrapper(self.calcs,fermi_surf=fermi_surface)
        ne=float(dos_range[1]-dos_range[0])/de

        calc_type='Density of States'
        print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                                            AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))

        if projected==True:
            AFLOWpi.scfuj.paopy_pdos_wrapper(self.calcs)
            calc_type='Projected Density of States'

            print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))

        if fermi_surface==True:
            calc_type='Fermi Surface'
            print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))

    def bands(self,nk=1000,band_topology=False,ipol=0,jpol=1,spol=2):
        '''
        Calculate energy eigenvalues along high symmetry path using PAOFLOW

        Keyword Arguments:
              nk (int): the approximate number of sampling points in the Brillouin Zone along
                         the entirety of the path between high symmetry points. Points are 
                         chosen so that they are equidistant along the entirety of the path.

              band_topology (bool): do the band topology calculation along the path
              ipol (int): first direction when calculating band topology quantites
              jpol (int): second direction when calculating band topology quantites
              spol (int): third direction when calculating band topology quantites

        Returns:
              None

        '''

        fermi_surface=False
        nbnd=None
        eShift=15.0
        cond_bands=True

        AFLOWpi.scfuj.paopy_bands_wrapper(self.calcs,band_topology=band_topology,fermi_surface=fermi_surface,ipol=ipol,jpol=jpol,spol=spol,nk=nk)

        calc_type='Electronic Band Structure'
        print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                                            AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))
        if band_topology:
            calc_type='Band Topology'
            print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))
        if fermi_surface:
            calc_type='Fermi Surface'
            print((AFLOWpi.run._colorize_message('ADDING TB STEP:  ',level='GREEN',show_level=False)+\
                AFLOWpi.run._colorize_message(calc_type,level='DEBUG',show_level=False)))


def _form_TB_dir(oneCalc,ID,from_ls=True):
    if from_ls:
        AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.save'])
#        AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.save/atom_proj.dat'],#%oneCalc['_AFLOWPI_PREFIX_']],
#                                         glob=True,first_node_only=True)
#        AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.save/data-file.xml'],#%oneCalc['_AFLOWPI_PREFIX_']],
#                                         glob=True,first_node_only=True)
    try:
        save_dir = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],oneCalc['_AFLOWPI_PREFIX_']+'.save')
        TB_dir = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_TB.save'%ID)
        if not os.path.exists(TB_dir):
            os.mkdir(TB_dir)
        data_file_dft = os.path.join(save_dir,'data-file.xml')
        atomic_proj_dat = os.path.join(save_dir,'atomic_proj.xml')
#        atomic_proj_dat = os.path.join(save_dir,'atomic_proj.dat')
        try:
            shutil.copy(data_file_dft,TB_dir)
        except:
            data_file_dft = os.path.join(save_dir,'data-file-schema.xml')
            shutil.copy(data_file_dft,TB_dir)

        shutil.copy(atomic_proj_dat,TB_dir)

        for UPF_FILE in glob.glob("%s/*.UPF"%(save_dir)):
            shutil.copy(UPF_FILE,TB_dir)
        for UPF_FILE in glob.glob("%s/*.upf"%(save_dir)):
            shutil.copy(UPF_FILE,TB_dir)


    except Exception as e:
        print(e)


def _run_want_bands(__submitNodeName__,oneCalc,ID,num_points=1000,cond_bands=True,compute_ham=False,proj_thr=0.95,proj_sh=5.5):
    nscf_ID=ID+'_nscf'

    Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
    eShift=float(Efermi)+10.0
    
    if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
        execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
    else:
        execPrefix=''

    want_dict = AFLOWpi.scfuj.WanT_bands(oneCalc,ID=ID,eShift=proj_sh,num_points=num_points,cond_bands=cond_bands,compute_ham=compute_ham,proj_thr=proj_thr)



    for want_ID,want_calc in list(want_dict.items()):
        AFLOWpi.run._oneRun(__submitNodeName__,want_calc,want_ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='custom',execPath='./want_bands.x',)

    AFLOWpi.prep._clean_want_bands(oneCalc,ID)

    return oneCalc,ID

def _rename_projectability(oneCalc,ID):
#    nspin = int(AFLOWpi.scfuj.chkSpinCalc(oneCalc,ID))
#    if nspin==2:
        proj_up = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'projectability_dn.txt')
        proj_dn = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'projectability_up.txt')
        try:
            os.rename(proj_dn,proj_dn_new)
        except:
            pass
        proj_dn_new = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_projectability_dn.txt'%ID)
        proj_up_new = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_projectability_up.txt'%ID)        
        try:
            os.rename(proj_up,proj_up_new)
        except:
            pass
        proj = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'projectability.txt')
        proj_new = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_projectability.txt'%ID)
        try:
            os.rename(proj,proj_new)
        except:
            pass

# def _run_want_eff_mass(__submitNodeName__,oneCalc,ID,temperature=[0,800],step=10):
#     nscf_ID=ID+'_nscf'

#     if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
#         execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
#     else:
#         execPrefix=''



#     if '__effmass_counter__' not in oneCalc.keys():
#         #if an old effective mass data file exists delete it before we start
#         effmass_datafile_by_temp = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_WanT_effmass.dat'%ID)
#         if os.path.exists(effmass_datafile_by_temp):
#             os.remove(effmass_datafile_by_temp)
#         #set counter to zero
#         oneCalc['__effmass_counter__']=0
#         AFLOWpi.prep._saveOneCalc(oneCalc,ID)

#     cell_params = AFLOWpi.retr._getCellParams(oneCalc,ID)
#     k_grid = AFLOWpi.prep.getMPGrid(cell_params,offset=True,string=False)
#     try:
#         k_grid = [int(float(x)*10.0) for x in k_grid.split()[:3]]
#     except:
#         k_grid=[20,20,20]

#     Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
#     eShift=float(Efermi)+10.0

#     step_holder=step
#     step = (float(temperature[1])-float(temperature[0])+float(step_holder))/float(step)
#     temps = numpy.linspace(float(temperature[0]),float(temperature[1]),step)

#     #some constants
#     h_bar = numpy.float64(1.05457180*10.0**-34.0)
#     k_b   = numpy.float64(1.38064852*10.0**-23.0)
#     m_e   = numpy.float64(9.10938356*10.0**-31.0)

#     sf = numpy.power(k_b*m_e/(2.0*numpy.pi*numpy.power(h_bar,2.0)),(3.0/2.0))#*numpy.power((1.0/cm2m),2.0)


#     for temp_step in range(len(temps)):
#         if temp_step<oneCalc['__effmass_counter__']:
#             continue

#         want_dos_calc = AFLOWpi.scfuj.WanT_dos(oneCalc,ID,k_grid=k_grid,pdos=False,boltzmann=False,eShift=eShift,cond_bands=True,temperature=temps[temp_step])
#         this_temp = '%8.4f ' % float(temps[temp_step])
#         for want_dos_ID,want_dos in want_dos_calc.iteritems():
#             AFLOWpi.run._oneRun(__submitNodeName__,want_dos,want_dos_ID,engine='espresso',calcType='custom',execPath='./effmass.x',execPrefix=execPrefix,execPostfix='')

#             effmass_datafile = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.dat'%want_dos_ID)
#             effmass_datafile_by_temp = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_WanT_effmass.dat'%ID)
#             with open(effmass_datafile,'r') as emdfo:
#                 data_by_line = emdfo.readlines()



#             with open(effmass_datafile_by_temp,'a+') as emdfo:

#                 for data_line in range(len(data_by_line)):
#                     if temp_step==0:
#                         if data_line==0:
#                             temp_as_str = 'Temperature '+data_by_line[data_line]

#                     else:
#                         if data_line==0:
#                             continue
#                         else:
#                             try:

#                                 emass = numpy.float64(data_by_line[data_line].strip('\n').split()[-1])


#                             except Exception as e:
#                                 print e
#                                 continue
#                             line_write = this_temp + data_by_line[data_line].strip('\n')+'\n'#+' %s\n'%(N_s)
#                             emdfo.write(line_write)

#         oneCalc['__effmass_counter__']=temp_step
#         AFLOWpi.prep._saveOneCalc(oneCalc,ID)

#     del oneCalc['__effmass_counter__']
#     AFLOWpi.prep._saveOneCalc(oneCalc,ID)

#     return oneCalc,ID

def _run_want_dos(__submitNodeName__,oneCalc,ID,dos_range=[-6,6],k_grid=None,project=True,num_e=2001,cond_bands=True,fermi_surface=False,compute_ham=False,proj_thr=0.95,proj_sh=5.5):
    nscf_ID=ID+'_nscf'

    if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
        execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
    else:
        execPrefix=''

    Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
#    eShift=float(Efermi)+10.0

    want_dos_calc = AFLOWpi.scfuj.WanT_dos(oneCalc,ID,energy_range=dos_range,k_grid=k_grid,pdos=project,boltzmann=False,num_e=num_e,eShift=proj_sh,cond_bands=cond_bands,fermi_surface=fermi_surface,compute_ham=compute_ham,proj_thr=proj_thr,)


    for want_dos_ID,want_dos in list(want_dos_calc.items()):
        AFLOWpi.run._oneRun(__submitNodeName__,want_dos,want_dos_ID,engine='espresso',calcType='custom',execPath='./want_dos.x',execPrefix=execPrefix,execPostfix='')

        if project==True:
            spin_state = want_dos_ID.split('_')[-1].strip()


    if len(list(want_dos_calc.keys()))>1:
        AFLOWpi.prep._combine_pol_pdos(oneCalc,ID)


    return oneCalc,ID


def _convert_tb_pdos(oneCalc,ID,spin=0):
    
        want_pdos_ext_glob = '_WanT_dos-*.dat'
        

        #change the TB pdos file names so that they can be read by sumpdos
        if spin == -1:
            spin_postfix='_down'
            dat_postfix ='_1'
        elif spin == 1:
            spin_postfix='_up'
            dat_postfix ='_0'
        else:
            spin_postfix=''
            dat_postfix ='_0'


        qe_pdos_out_str = AFLOWpi.retr._getOutputString(oneCalc,ID+'_pdos')
        rename_info_re = re.compile(r'state #\s*(\d*): atom\s*(\d+)\s*\(\s*(\S*)\s*\).*wfc\s*(\d+)\s*\(l=(\d+).*\)\n')

        #first check the QE projwfc.x output for the orbital
        #and species label for each state #

        state_info_list = rename_info_re.findall(qe_pdos_out_str)
        #if it found the info on the states by their numbers
        if len(state_info_list)==0:
            rename_info_re = re.compile(r'state #\s*(\d*): atom\s*(\d+)\s*\(\s*(\S*)\s*\).*wfc\s*(\d+).*l=(\d+).*m_j=([\s-][.\d]+).*\n')
            state_info_list = rename_info_re.findall(qe_pdos_out_str)

        pdos_dict=OrderedDict()
        if len(state_info_list)!=0:
            for i in range(len(state_info_list)):

                    state_num = int(state_info_list[i][0].strip())
                    atom_num  = state_info_list[i][1].strip()
                    atom_spec = state_info_list[i][2].strip()
                    wfc_num   = state_info_list[i][3].strip()
                    orb_l     = int(state_info_list[i][4].strip())
                    try: 
                        orb_m_j     = state_info_list[i][5].strip().replace('-','m')
                    except: orb_m_j=-999
                    #translate atomic number "l" to orbital name (i.e. s,p,d,f)
                    orb_type=['s','p','d','f']

                    #the orig file name from WanT output
                    orig_name = '%d_pdosdk%s.dat'%(state_num-1,dat_postfix)

                    orig_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],orig_name)
                    if not os.path.exists(orig_path):
                        orig_name = '%d_pdos%s.dat'%(state_num-1,dat_postfix)
                        orig_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],orig_name)
                    #the renamed filename
                    if orb_m_j==-999:
                        new_name = '%s_TB%s.pdos_atm#%s(%s)_wfc#%s(%s)'%(ID,spin_postfix,
                                                                         atom_num,atom_spec,
                                                                         wfc_num,orb_type[orb_l])
                    else:
                        new_name = '%s_TB%s.pdos_atm#%s(%s)_wfc#%s(%s)_%s'%(ID,spin_postfix,atom_num,
                                                                            atom_spec,wfc_num,
                                                                            orb_type[orb_l],orb_m_j)

                    
                    try:
                        pdos_dict[new_name].append(orig_name)
                    except Exception as e:
                        pdos_dict[new_name]=[orig_name]

        for k,v in list(pdos_dict.items()):

            old_name_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],v[0])
            dat = numpy.loadtxt(old_name_path)                    
            for i in range(1,len(v)):
                old_name_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],v[i])
                dat[:,1] += numpy.loadtxt(old_name_path)[:,1]
            new_name_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],k)
            numpy.savetxt(new_name_path,dat)






def _combine_pol_pdos(oneCalc,ID):
    glob_ID =  AFLOWpi.prep._return_ID(oneCalc,ID,step_type='PAO-TB',last=True,straight=False)
    glob_ID +='_TB'

    glob_ID_up=glob_ID+'_up'
    glob_ID_dn=glob_ID+'_down'

    subdir=oneCalc['_AFLOWPI_FOLDER_']

    pdos_files_up = sorted(glob.glob(os.path.join(subdir,'%s.pdos_atm*' % (glob_ID_up))))
    pdos_files_dn = sorted(glob.glob(os.path.join(subdir,'%s.pdos_atm*' % (glob_ID_dn))))
                             
    for pdos_file in range(len(pdos_files_up)):
        output_list=[]
        pdos_file_up = pdos_files_up[pdos_file]
        pdos_file_dn = pdos_files_dn[pdos_file]

        with open(pdos_file_up) as pdfuo:
            pdos_string_up = pdfuo.readlines()
        with open(pdos_file_dn) as pdfdo:
            pdos_string_dn = pdfdo.readlines()

        for entry in range(len(pdos_string_up)):
            try:
                entry_up = pdos_string_up[entry].split()
                entry_dn = pdos_string_dn[entry].split()
                energy = entry_up[0]
                val_up = entry_up[1]
                val_dn = entry_dn[1]
            
                output_list.append('%s %s %s' % (energy,val_up,val_dn))
            except:
                pass



        output_str = '\n'.join(output_list)

        input_file_name = os.path.basename(pdos_file_up)
        input_file_name_split = input_file_name.split('.')[-1]
        output_file_name = ID+'_TB.'+input_file_name_split
        output_file_path = os.path.join(subdir,output_file_name)        

        with open(output_file_path,'w') as pdfco:
            pdfco.write(output_str)


class tb_plotter:
        '''
        Class for adding common plotting functions from AFLOWpi.plot module to the high level user 
        interface. 

        '''
        def __init__(self,calcs):
                self.calcs=calcs

        def topology(self,en_range=[-5,5],runlocal=False,postfix=''):
            '''
            Plot band topology data generated by PAOFLOW

            Keyword Arguments:
                  en_range (list): a tuple or list of the range of energy around the fermi/Highest
                            occupied level energy that is to be included in the plot.
                  runlocal (bool): if True, run plotting routine from user script. 
                                   Useful if replotting previously generated data with different en_range
                  postfix (str): a string of an optional postfix to the plot filename for every
                               calculation.

            Returns:
                  None

            '''

            AFLOWpi.plot.band_topology(self.calcs,yLim=en_range,DOSPlot='',runlocal=runlocal,postfix=postfix,tight_banding=False)

            calc_type='Plot Band Topology'
            print(('                 %s'% (calc_type)))


        def ahc(self,en_range=None,runlocal=False):
            '''
            Plot Anomalous Hall conductivity generated by PAOFLOW

            Keyword Arguments:
                  en_range (list): energy range around fermi level for plot. Default is the range of energy in which AHC was calculated
                  runlocal (bool): if True, run plotting routine from user script. 
                                   Useful if replotting previously generated data with different en_range
            Returns:
                  None

            '''


            if runlocal:
                try:
                    AFLOWpi.plot.__plot_berry_cond(oneCalc,ID,spin=False,en_range=en_range)
                except: pass
                try:
                    AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=False,real=False,en_range=en_range)
                except: pass
                try:
                    AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=False,real=True,en_range=en_range)
                except: pass

            addit = "AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=False,real=False,en_range=%s)"%en_range
            AFLOWpi.prep.addToAll_(self.calcs,"PLOT",addit)
            addit = "AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=False,real=True,en_range=%s)"%en_range
            AFLOWpi.prep.addToAll_(self.calcs,"PLOT",addit)
            addit = "AFLOWpi.plot.__plot_berry_cond(oneCalc,ID,spin=False,en_range=%s)"%en_range
            AFLOWpi.prep.addToAll_(self.calcs,"PLOT",addit)

            calc_type='Plot Anomalous Hall Conductivity'
            print(('                 %s'% (calc_type)))

        def shc(self,en_range=None,runlocal=False):
            '''
            Plot spin Hall conductivity generated by PAOFLOW

            Keyword Arguments:
                  en_range (list): energy range around fermi level for plot. Default is the range of energy in which SHC was calculated
                  runlocal (bool): if True, run plotting routine from user script. 
                                   Useful if replotting previously generated data with different en_range
            Returns:
                  None

            '''
            if runlocal:
                try:
                    AFLOWpi.plot.__plot_berry_cond(oneCalc,ID,spin=True,en_range=en_range)
                except: pass
                try:
                    AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=True,real=False,en_range=en_range)
                except: pass
                try:
                    AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=True,real=True,en_range=en_range)
                except: pass



            addit = "AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=True,real=False,en_range=%s)"%en_range
            AFLOWpi.prep.addToAll_(self.calcs,"PLOT",addit)
            addit = "AFLOWpi.plot.__plot_dichroism(oneCalc,ID,spin=True,real=True,en_range=%s)"%en_range
            AFLOWpi.prep.addToAll_(self.calcs,"PLOT",addit)
            addit = "AFLOWpi.plot.__plot_berry_cond(oneCalc,ID,spin=True,en_range=%s)"%en_range
            AFLOWpi.prep.addToAll_(self.calcs,"PLOT",addit)

            calc_type='Plot Spin Hall Conductivity'
            print(('                 %s'% (calc_type)))

        def opdos(self,en_range=[-5,5],runlocal=False,postfix=''):
            '''
            Plot orbital projected DOS generated by PAOFLOW

            Keyword Arguments:
                  en_range (list): list of upper and lower bounds around fermi energy for plot range
                  runlocal (bool): if True, run plotting routine from user script. 
                                   Useful if replotting previously generated data with different en_range
                  postfix (str): postfix to plot file names
            Returns:
                  None

            '''

            AFLOWpi.plot.opdos(self.calcs,yLim=en_range,runlocal=runlocal,postfix=postfix,tight_binding=True)

            calc_type='Plot Orbital Projected DOS of PAO-TB Representation'
            print(('                 %s'% (calc_type)))


        def apdos(self,en_range=[-5,5],runlocal=False,postfix=''):
            '''
            Plot atom projected DOS generated by PAOFLOW

            Keyword Arguments:
                  en_range (list): list of upper and lower bounds around fermi energy for plot range
                  runlocal (bool): if True, run plotting routine from user script. 
                                   Useful if replotting previously generated data with different en_range
                  postfix (str): postfix to plot file names
            Returns:
                  None

            '''

            AFLOWpi.plot.apdos(self.calcs,yLim=en_range,runlocal=runlocal,postfix=postfix,tight_binding=True)

            calc_type='Plot Atom Projected DOS of PAO-TB Representation'
            print(('                 %s'% (calc_type)))


        def transport(self,runlocal=False,postfix='',en_range=None):
                '''
                Plot Boltzmann transport properties calculated with PAOFLOW

                Keyword Arguments:
                      en_range (list): list of upper and lower bounds around fermi energy for plot range
                      runlocal (bool): if True, run plotting routine from user script. 
                                       Useful if replotting previously generated data with different en_range
                      postfix (str): postfix to plot file names                
                Returns:
                      None

                '''

                AFLOWpi.plot.transport_plots(self.calcs,runlocal=runlocal,postfix=postfix,x_range=en_range)
                
                calc_type='Plot Boltzmann Transport'
                print(('                 %s'% (calc_type)))


        def optical(self,runlocal=False,postfix='',en_range=None):
                '''
                Plot optical properties calculated with PAOFLOW

                Keyword Arguments:
                      en_range (list): list of upper and lower bounds around fermi energy for plot range
                      runlocal (bool): if True, run plotting routine from user script. 
                                       Useful if replotting previously generated data with different en_range
                      postfix (str): postfix to plot file names                
                Returns:
                      None

                '''



                AFLOWpi.plot.optical_plots(self.calcs,runlocal=runlocal,postfix=postfix,x_range=en_range)
                
                calc_type='Plot Optical Epsilon'
                print(('                 %s'% (calc_type)))



        def bands(self,en_range=[-5,5],DOSPlot='',runlocal=False,postfix=''):
            '''
            Plot electronic band structure generated by calculation engine

            Keyword Arguments:
                  en_range (list): a tuple or list of the range of energy around the fermi/Highest
                            occupied level energy that is to be included in the plot.
                  DOSPlot (str): a string that flags for the option to have either a DOS plot
                               share the Y-axis of the band structure plot. Options include:

                               ""      | A blank string will cause No Density of
                                       | States plotted alongside the Band Structure

                               "APDOS" | Atom Projected Density of States
                               "DOS"   | Normal Density of States


                  runlocal (bool): if True, run plotting routine from user script. 
                                   Useful if replotting previously generated data with different en_range

                  postfix (str): a string of an optional postfix to the plot filename for every
                               calculation.

            Returns:
                  None

            '''

            AFLOWpi.plot.bands(self.calcs,yLim=en_range,DOSPlot=DOSPlot,runlocal=runlocal,postfix=postfix,tight_banding=True)

            calc_type='Plot Electronic Band Structure of PAO-TB Representation'
            if DOSPlot=='DOS':
                    calc_type+=' with Density of States'
            if DOSPlot=='APDOS':
                    calc_type+=' with APDOS'
            print(('                 %s'% (calc_type)))

        # def dos(self,en_range=[-5,5],runlocal=False,postfix=''):
        #     pass


def _run_tb_ham_prep(__submitNodeName__,oneCalc,ID,config=None,kp_factor=2.0,cond=1,ovp=False,band_factor=1.25,tetra_nscf=False,wsyminv=False):
        execPrefix = ''
        execPostfix = ''
        oneCalcID = ID

        try:
            pw_path = os.path.join(AFLOWpi.prep._ConfigSectionMap('prep','engine_dir'),'pw.x')
            shutil.copy(pw_path,oneCalc['_AFLOWPI_FOLDER_'])
        except: pass

        def abortIFRuntimeError(subdir, ID):
            with open(os.path.join(subdir, "%s.out"%ID)) as ifo:
                outfile =  ifo.read()
 
            errorList = re.findall(r'from (.*) : error #.*\n',outfile)
            if len(errorList) > 0:        
                logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
                print(("Error in %s.out -- ABORTING ACBN0 LOOP"%ID))                    
                raise SystemExit



        if '__runList__' not in list(oneCalc.keys()):
            oneCalc['__runList__']=[]

            
        if config is not None:
                AFLOWpi.prep._forceGlobalConfigFile(config)
                logging.debug('forced config %s' % config)
        else:
                try:
                        config = AFLOWpi.prep._getConfigFile()
                        AFLOWpi.prep._forceGlobalConfigFile(config)
                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)


        if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
            execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

        else:
            execPrefix=''


        if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
                execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
        else:
                execPostfix=''


        if AFLOWpi.prep._ConfigSectionMap('run','engine') == '':
                engine = AFLOWpi.prep._ConfigSectionMap('run','engine')
        else:
                engine = 'espresso'


        subdir = oneCalc['_AFLOWPI_FOLDER_']
        oneCalc['_AFLOWPI_CONFIG_']=config

        if 'scf' not in oneCalc['__runList__']:

            try:
                npool=AFLOWpi.retr._get_pool_num(oneCalc,ID)        

                if npool!=1:
                    if len(re.findall(r'npool[s]*\s*(?:\d*)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool[s]*\s*(?:\d*)','npool %s'%npool,execPostfix)
                        logging.debug(execPostfixPrime)

            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)


##################################################################################################################
            AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)


            oneCalc['__runList__'].append('scf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            



            nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_factor,unoccupied_states=cond,band_factor=band_factor,tetra_nscf=tetra_nscf,wsyminv=wsyminv)  



        else:
            '''if we are restarting from a job killed from going walltime 
            try to load ID_nscf and if we can't then just make a new one'''
            try:
                nscf_ID='%s_nscf' % ID
                nscf_calc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],nscf_ID)                
                '''we have to make sure nscf step has the correct walltime and start time if it's a restart'''
                nscf_calc['__walltime_dict__']=oneCalc['__walltime_dict__']
            except Exception as e:
                try:
                    nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_factor,
                                                                      band_factor=band_factor)  

                except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)

        nscf_exec_postfix = execPostfix_LOCAL = AFLOWpi.prep._ConfigSectionMap('TB','exec_postfix_nscf')            
        if nscf_exec_postfix != "":
            execPostfix = nscf_exec_postfix 
        
##################################################################################################################
        if 'nscf' not in oneCalc['__runList__']:


            try:
                npool=AFLOWpi.retr._get_pool_num(nscf_calc,nscf_ID)        

                if npool!=1:
                    if len(re.findall(r'npool\s*(?:\d+)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool\s*(?:\d+)','npool %s'%npool,execPostfix)
                        logging.debug(execPostfixPrime)

            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

            AFLOWpi.run._oneRun(__submitNodeName__,nscf_calc,nscf_ID,execPrefix=execPrefix,
                                execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)
            AFLOWpi.retr._writeEfermi(nscf_calc,nscf_ID)

            abortIFRuntimeError(subdir, nscf_ID)

            oneCalc['__runList__'].append('nscf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)       
##################################################################################################################
        pdos_calc,pdos_ID = AFLOWpi.scfuj.projwfc(oneCalc,ID,paw=False,ovp=ovp)


        if not re.match('northo',execPostfix) or not re.match('no',execPostfix):
            execPostfix+=' -northo 1'

        execPrefix_LOCAL = AFLOWpi.prep._ConfigSectionMap('run','exec_prefix')
        execPostfix_LOCAL = AFLOWpi.prep._ConfigSectionMap('run','exec_postfix')            
        
        if 'pdos' not in oneCalc['__runList__']:
            pdosPath = os.path.join(AFLOWpi.prep._ConfigSectionMap('prep','engine_dir'),'projwfc.x')

            shutil.copy(pdosPath,oneCalc['_AFLOWPI_FOLDER_'])

            pdosPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'projwfc.x')

            AFLOWpi.run._oneRun(__submitNodeName__,pdos_calc,pdos_ID,execPrefix=execPrefix,
                                execPostfix=execPostfix,engine='espresso',calcType='custom',
                                executable='projwfc.x',execPath=pdosPath)
#############
            oneCalc['__runList__'].append('pdos')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            abortIFRuntimeError(subdir, pdos_ID)



        eFermi=0.0


        AFLOWpi.prep._form_TB_dir(oneCalc,ID)
        eFermi=10.0

        splitInput = AFLOWpi.retr._splitInput(nscf_calc['_AFLOWPI_INPUT_'])
        del oneCalc['__runList__']

        dos_fermi = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_PAOFLOW_dos.efermi'%ID)

        with open(dos_fermi,'w') as ifo:
                ifo.write(str(0.0))

        return oneCalc,ID




def _clean_want_bands(oneCalc,ID):


    try:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_bands.out'%ID)[-1]
    except:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_bands_up.out'%ID)[-1]

    with open(want_stdout_path,'r') as in_file_obj:
        in_string = in_file_obj.read()

    path = AFLOWpi.retr._getPath(0.01,oneCalc,ID=ID)

    plot_bool=[]
    path_name=[]
    path_split = [x for x in  path.split('\n')[1:] if len(x.strip())]
    for i in path_split:
        path_name.append(i.split()[-1])
        if  int(i.split()[3]):

            plot_bool.append(True)
        else:
            plot_bool.append(False)

    gg = re.findall("  Number of kpts in each segment\n((?:.*:\W+(?:\d*)\n)*)",in_string)
#    num = [int(x) for x in re.findall('line.*:\W+(\d+)',gg[0])]
    num = [int(x) for x in re.findall('line\s*\d+:\s*(\d+)\s*\n',in_string)]
    total = 0
    include=[]
    #print path_name
    output_path_string = ''
    for i in range(len(num)):
        total+=num[i]+1
        try:
            if plot_bool[i]:
                if i==0:
                    output_path_string+='%s %s\n' %(path_name[i],num[i])
                else:
                    output_path_string+='%s %s\n' %(path_name[i],num[i]+1)
            else:
                output_path_string+='%s %s\n' %(path_name[i],0)
            for j in range(num[i]):
                include.append(plot_bool[i])

            include.append(True)
        except:
            pass
    print(include)
    #    print print_out
    output_path_string+='%s %s' %(path_name[-1],0)+'\n' 
    
    want_bands_data_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_bands_want.dat'%ID)
    if os.path.exists(want_bands_data_path):
        data_file_list=[want_bands_data_path]
    else:
        want_bands_data_path_up = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_bands_want_up.dat'%ID)
        want_bands_data_path_dn = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_bands_want_down.dat'%ID)
        data_file_list=[want_bands_data_path_up,want_bands_data_path_dn,]



    ret_data=[]
    for want_bands_data_path in data_file_list:
        with open(want_bands_data_path,'r') as in_file_obj:
            bands_dat = in_file_obj.read()
        split_bands = bands_dat.split('\n')

        split_data = []
        per_band=[]
        for i in split_bands:
        #    print len(i.strip())
            if not len(i.strip()):
                if len(per_band)!=0:
                    split_data.append(per_band)
                per_band=[]
            else:
                per_band.append(i)


        final_data=''
        for i in range(len(split_data)):
            for j in range(len(split_data[i])):
                if include[j]:
                    final_data+= split_data[i][j]+'\n'

            final_data+='\n'
            
        ret_data.append(final_data)
        want_bands_data_path_new = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],want_bands_data_path[:-4]+'_cleaned.dat')
        print(want_bands_data_path_new)
#        want_bands_data_path_new=want_bands_data_path
        with open(want_bands_data_path_new,'w') as in_file_obj:
            in_file_obj.write(final_data)

#    return ret_data
#def _get_pp_nawf(oneCalc,ID):
    
