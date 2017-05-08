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
import copy 

def _collect_fd_field_forces(oneCalc,ID,for_type='raman'):
    """
    Runs at the end of the scf calculations for the Finite Difference phonon
    calculations. It uses regex to pull the forces to be saved to files for
    parsing by fd_ifc.x
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
         for_type (str): type of file to pull from in FD_PHONON dir (choices are 'raman', 'born', 'epol')

    Returns:
         force_out_str (str) a string of the contents of the file
          
    """

    if for_type=='raman':
        val_type='force'
        num=range(19)
    elif for_type=='born':
        val_type='force'
        num=[0,2,4,6,1,3,5]
    elif for_type=='epol':
        val_type='epol'
        num=[0,2,4,6,1,3,5]

    force_out_str=''
    for index in num:
        try:
            force_file_path=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'./FD_PHONON/%s.%02d'%(val_type,index))
            with open(force_file_path,'r') as force_file:
                force_str=force_file.read()
                force_out_str+=force_str
        except:
            pass

    return force_out_str

def _pull_polarization(oneCalc,ID):
    """
    Runs at the end of the scf calculations for the Finite Difference phonon
    calculations. It uses regex to pull the forces to be saved to files for
    parsing by fd_ifc.x
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
         None

    Returns:
         None 
          
    """

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out'%ID),'r') as out_file_obj:
        out_string = out_file_obj.read()
    
    e_pol_regex = re.compile('Electronic Dipole on Cartesian axes\s*\n((?:\W*[123]\s*.*\n)+)')
    e_ion_regex = re.compile('Ionic Dipole on Cartesian axes\s*\n((?:\W*[123]\s*.*\n)+)')

    try:
        ele_block=e_pol_regex.findall(out_string)[-1]
        ion_block=e_ion_regex.findall(out_string)[-1]



    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        return
    e_pol_split = [float(x.split()[-1]) for x in ele_block.split('\n') if len(x.strip())!=0]
    e_ion_split = [float(x.split()[-1]) for x in ion_block.split('\n') if len(x.strip())!=0]
    
    
    e_tot_split = [e_pol_split[0],e_pol_split[1],e_pol_split[2],]

    epol_out_string=''

    epol_out_string+='%25.15e%25.15e%25.15e\n'%(float(e_tot_split[0]),float(e_tot_split[1]),float(e_tot_split[2]),)

    force_postfix='.'.join(ID.split('.')[1:])
    force_postfix=force_postfix.split('_')[0]

    with open(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'epol.%s'%force_postfix),'w+') as out_file_obj:
        out_file_obj.write(epol_out_string)


def _pull_eps_out(oneCalc,ID):
    """
    Runs at the end of the scf calculations for the Finite Difference phonon
    calculations. It uses regex to pull the eps_0 to be saved to files for
    use by matdyn.x
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
         None

    Returns:
         eps_string (str): string of the eps
          
    """

    eps_string=''
    eps_regex=re.compile('Dielectric tensor\n.*\n\s*([0-9\s.-]*)\n')
    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_epol.out'%ID),'r') as out_file_obj:
            out_string = out_file_obj.read()

        eps_string=eps_regex.findall(out_string)[-1]
        
        
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    return eps_string

def _pull_born_out(oneCalc,ID):
    """
    Runs at the end of the scf calculations for the Finite Difference phonon
    calculations. It uses regex to pull the born effective charges to be saved 
    to files for use by matdyn.x
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
         None

    Returns:
         born_string (str) string of the born eff. charges
          
    """

    born_string=''
    born_regex=re.compile('atom\s*\d*\n(?:([-.0-9\s]+)+)')
    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_born.out'%ID),'r') as out_file_obj:
            out_string = out_file_obj.read()
        born_charges=born_regex.findall(out_string)


        for i in range(len(born_charges)):
            num_index=i+1
            born_string+='         %s\n'%num_index
            born_string+=born_charges[i]

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    return born_string

def _gen_fd_input(oneCalc,ID,for_type='raman',de=0.003):
    """
    Generates the input files for the finite field calcs from the input file string in oneCalc
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
         for_type (str): which type of FD input do we want to generate files for. 'raman','born','eps'
         de (float): applied electric field strength

    Returns:
         None
          
    """

    if for_type=='raman':
        card='RAMAN_TENSOR'
        input_name_type='raman'
        npol=4
    elif for_type=='born':
        card='BORN_CHARGES'        
        input_name_type='zeu'
        npol=2
    elif for_type=='epol':
        card='DIELECTRIC_TENSOR'
        input_name_type='eps'
        npol=2


    header='''&inputfd
   prefix='%s'
   de_%s=%s
   npol_%s=%s
/
'''% (oneCalc['_AFLOWPI_PREFIX_'],input_name_type,de,input_name_type,npol)



    header+=card+'\n'

    forces = AFLOWpi.run._collect_fd_field_forces(oneCalc,ID,for_type=for_type)
    header+=forces
    
    finite_field_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_%s.in'%(ID,for_type))
    with open(finite_field_input,'w') as force_file:
        force_file.write(header)

def _swap_scf_inputs(oneCalc,ID):
    """
    Swaps the value of oneCalc['_AFLOWPI_INPUT_'] in oneCalc to the finite field input
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          oneCalc (dict): dictionary of one of the calculations
          
    """

    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'r') as in_by_ID_obj:
            input_string = in_by_ID_obj.read()
        oneCalc['_AFLOWPI_INPUT_']=input_string
    except:
        pass

    return oneCalc

def _setup_raman(calc_subset,orig_oneCalc,orig_ID,field_strength=0.001,field_cycles=3,for_type='raman'):
    """
    Sets up the input files and the command to run the finite field calculations.
    
    Arguments:
          calc_subset (dict): dictionary of dictionaries of the FD_PHONON calculations 
          orig_oneCalc (dict): oneCalc of one calculation the main calc set

    Keyword Arguments:
         feild_strength (float): applied electric field strength
         field_cycles (int) number of berry phase cycles
         for_type (str): which type of FD input do we want to generate files for. 'raman','born','eps'

    Returns:
         None
          
    """

    input_string=orig_oneCalc['_AFLOWPI_INPUT_']
    subset_keys = calc_subset.keys()

    in_dict=AFLOWpi.retr._splitInput(input_string)

    #if possibly a metal don't do LOTO

    if 'degauss' not in  in_dict['&system'].keys():
        pass
    else:
        bandgap_type = AFLOWpi.retr.gap_type(orig_oneCalc,orig_ID)
        if bandgap_type not in ['p-type','n-type','insulator']:
            return
        else:
            if 'degauss' in in_dict['&system'].keys():
                del in_dict['&system']['degauss']
            else:
                pass

            if "occupations" in in_dict['&system'].keys():
                val=eval(in_dict['&system']['occupations'].lower())
                if val=="smearing":
                    del in_dict['&system']['occupations']
                else:
                    pass
            else:
                pass

        input_string=AFLOWpi.retr._joinInput(in_dict)


    #copy for when we append finite field calcs to the calclog
    calc_subset_extended=copy.deepcopy(calc_subset)

    if for_type=='raman':
        val_type='force'
        num=range(19)
    elif for_type=='born':
        val_type='force'
        num=[0,2,4,6,1,3,5]
    elif for_type=='epol':
        val_type='epol'
        num=[0,2,4,6,1,3,5]

    num_calcs=19
    subset_size=len(calc_subset)
    
    for index in range(len(num)):        
        fd_input = AFLOWpi.run._field_factory(input_string,field_strength=field_strength,for_type=for_type,field_cycles=field_cycles,dir_index=num[index])

        subset_key_index=num[index]%subset_size
        exec_counter_index = (index/subset_size)+1
        ID = subset_keys[subset_key_index]
        oneCalc=calc_subset[ID]

        #set prefix in finite field input to that of the FD phonon calc's input
        split_input = AFLOWpi.retr._splitInput(fd_input)
        split_input['&control']['prefix']='"%s"'%calc_subset[ID]['_AFLOWPI_PREFIX_']
        fd_input = AFLOWpi.retr._joinInput(split_input)

        FD_ID='finite_field.%02d'%num[index]
        file_path_fd_field_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%FD_ID)
    

        with open(file_path_fd_field_in,'w') as fd_in_obj:
            fd_in_obj.write(fd_input)


        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"""oneCalc_ff = AFLOWpi.run._swap_scf_inputs(oneCalc,'%s')""" % FD_ID)       
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','''AFLOWpi.prep._to_local_scratch(oneCalc_ff,'%s')'''%FD_ID)


        execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
        run_command="""if oneCalc['__execCounter__']<=%s:
        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s',execPrefix='%s',execPostfix='-northo 1',engine='espresso',calcType='scf',exit_on_error=True)
        oneCalc['__execCounter__']+=1
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
"""%(str(exec_counter_index),FD_ID,execPrefix)

        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',run_command)       

#        AFLOWpi.run._onePrep(oneCalc,ID,engine='espresso',calcType='scf',execPrefix=execPrefix,execPostfix=' -northo 1 ',alt_ID=FD_ID,)


        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"""AFLOWpi.run._pull_forces(oneCalc_ff,'%s')""" % FD_ID)       
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"""AFLOWpi.run._pull_polarization(oneCalc_ff,'%s')""" % FD_ID)       


        # add finite field to calcs for calclog
        oneCalc_FD=copy.deepcopy(oneCalc)
        oneCalc_FD['_AFLOWPI_INPUT_']=fd_input
        calc_subset_extended[FD_ID]=oneCalc_FD

    #update the calclogs to include FD_fields so when all phonon supercell
    #are finished but not finite fields the main script doesn't try to run
    chain_index=orig_oneCalc['__chain_index__']

#FIX THIS
#FIX THIS
    chain_logname='step_01'
###    chain_logname='step_%02d'% chain_index
#FIX THIS
#FIX THIS

#    AFLOWpi.prep.updatelogs(calc_subset_extended,chain_logname,runlocal=True)

def _field_factory(input_str,field_strength=0.001,for_type='raman',field_cycles=3,dir_index=0):
    """
    Modifies the QE pwscf input for the finite field calc of a given index
    
    Arguments:
         input_str (str): string of QE pwscf input file

    Keyword Arguments:
         feild_strength (float): applied electric field strength
         field_cycles (int) number of berry phase cycles
         for_type (str): which type of FD input do we want to generate files for. 'raman','born','eps'
         dir_index (int): index of the finite field direction to choose

    Returns:
         new_input (str): the modified input string
          
    """

    if for_type=='eps' or for_type=='born':
        field_dir=[[ 0.0, 0.0, 0.0,],  # 0
                   [-1.0, 0.0, 0.0,],  # -x
                   [ 1.0, 0.0, 0.0,],  #  x
                   [ 0.0,-1.0, 0.0,],  # -y
                   [ 0.0, 1.0, 0.0,],  #  y
                   [ 0.0, 0.0,-1.0,],  # -z
                   [ 0.0, 0.0, 1.0,],]  #  z
    if for_type=='raman':
        field_dir=[[ 0.0, 0.0, 0.0,],  # 0
                   [-1.0, 0.0, 0.0,],  # -x
                   [ 1.0, 0.0, 0.0,],  #  x
                   [ 0.0,-1.0, 0.0,],  # -y
                   [ 0.0, 1.0, 0.0,],  #  y
                   [ 0.0, 0.0,-1.0,],  # -z
                   [ 0.0, 0.0, 1.0,],  #  z
                   [-1.0,-1.0, 0.0,],  # -x-y
                   [ 1.0, 1.0, 0.0,],  #  xy
                   [ 1.0,-1.0, 0.0,],  #  x-y 
                   [-1.0, 1.0, 0.0,],  # -xy
                   [-1.0, 0.0,-1.0,],  # -x-z
                   [ 1.0, 0.0, 1.0,],  #  xz
                   [ 1.0, 0.0,-1.0,],  #  x-z
                   [-1.0, 0.0, 1.0,],  # -xz
                   [ 0.0,-1.0,-1.0,],  # -y-z
                   [ 0.0, 1.0, 1.0,],  #  yz
                   [ 0.0, 1.0,-1.0,],  #  y-z
                   [ 0.0,-1.0, 1.0,],]  #  -yz


    in_dict=AFLOWpi.retr._splitInput(input_str)
    in_dict['&control']['lelfield']='.TRUE.'
    in_dict['&control']['calculation']="'scf'"
    in_dict['&control']['nberrycyc']=field_cycles
    in_dict['&control']['tprnfor']='.TRUE.'
    in_dict['&control']['prefix']="'_finite_field.%02d'"%dir_index

    input_list=[]

#    for i in range(len(field_dir)):
    for j in range(len(field_dir[dir_index])):
        in_dict['&electrons']['efield_cart(%s)'%(j+1)]=field_dir[dir_index][j]*field_strength

    new_input = AFLOWpi.retr._joinInput(in_dict)

        
    return new_input
