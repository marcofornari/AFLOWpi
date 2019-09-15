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
import re
import os
import logging
import shutil
import pickle
import glob
import fnmatch
import copy 
import contextlib
import sys
import io
import subprocess


def _grab__linux_version(oneCalc,ID):
    status = subprocess.call("lsb_release -a", shell=True)
    print(status)
        
def _outputRegex(oneCalc,ID,regex=''):
        outFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)
        try:
                with open(outFileName,'r') as outFileObj:
                        outFileString = outFileObj.read()
                outRegex = re.compile(regex)
        except Exception as e:

                return ''
        try:
                return outRegex.findall(outFileString)[-1]

        except Exception as e:

                return ''

def _grab__auid(oneCalc,ID):
        outFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)
        inFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % ID)
        hashText=''
        try:
                with open(inFileName,'r') as inFileObj:
                        inFileText=inFileObj.read()

                with open(outFileName,'r') as outFileObj:
                        outFileText=outFileObj.read()

                auidHash = AFLOWpi.prep._hash64String(inFileText+outFileText)
                return auidHash
        except:
                return ''

def _grab__code(oneCalc,ID):
        outString = ''
        codeRegex=r'Program PWSCF ([A-Za-z0-9.]*)\s*'
        outString = AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=codeRegex)
        outString='qe'+outString[1:]
        return outString

def _grab__Bravais_lattice_orig(oneCalc,ID):
        outString = ''
        ibravRegex = r'bravais-lattice index\s*=\s*([-0-9]*)'
        ibrav = int(AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=ibravRegex))
        alat, cellOld = AFLOWpi.retr._getCellParams(oneCalc,ID) 

                    
        if(int(ibrav)==1):
                outString='CUB'

        elif(int(ibrav)==2):
                outString='FCC'

        elif(int(ibrav)==3):
               outString='BCC'

        elif(int(ibrav)==4):
                outString='HEX'

        elif(int(ibrav)==5):
                tx = cellOld[0][0]
                c = (((tx**2)*2)-1.0)
                c = -c
                alpha = numpy.arccos(c)

                if alpha < numpy.pi/2.0:
                    outString='RHL1'

                elif alpha > numpy.pi/2.0:
                    outString='RHL2'

        elif(int(ibrav)==6):
               outString='TET'

        elif(int(ibrav)==7):

                a = cellOld[1][0]*2
                c = cellOld[1][2]*2

                if(c < a):
                        outString='BCT1'

                elif(c > a):
                        outString='BCT2'

        elif(int(ibrav)==8):
                    outString='ORC'

        elif(int(ibrav)==9):
               outString='ORCC'

        elif(int(ibrav)==10):
                a1 = cellOld[0][0]*2
                b1 = cellOld[1][1]*2
                c1 = cellOld[2][2]*2
                myList = [a1, b1, c1]
                c = max(myList)
                a = min(myList)
                myList.remove(a)
                myList.remove(c)
                b = myList[0]

                if(1.0/a**2 > 1.0/b**2 + 1.0/c**2):
                        outString='ORCF1'

                elif(1.0/a**2 < 1.0/b**2 + 1.0/c**2):
                        outString='ORCF2'

                elif(1.0/a**2 == 1.0/b**2 + 1.0/c**2):
                        outString='ORCF3'

        elif(int(ibrav)==11):
                    outString='ORCI'

        elif(int(ibrav)==12):
                    outString='MCL'

        '''
        FINISH THIS
        FINISH THIS
        FINISH THIS
        FINISH THIS
        FINISH THIS
        FINISH THIS
        FINISH THIS
        FINISH THIS
        '''
                    
        return outString

def _grab__compound(oneCalc,ID):
        outString = ''
        outString = AFLOWpi.retr._getStoicName(oneCalc)

        return outString

def _grab__density(oneCalc,ID):
        outString = ''

        return outString


def _grab__eentropy_cell(oneCalc,ID):
        outString = ''

        return outString

def _grab__Egap(oneCalc,ID):
        outString = ''

        return outString

def _grab__energy_cell(oneCalc,ID):
        return  AFLOWpi.retr.grabEnergy(oneCalc,ID)

def _grab__position_labels(oneCalc,ID):
        return ','.join(AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_']))

def _grab__energy_cutoff(oneCalc,ID):
        outString = ''
        energyCutoffRegex = r'kinetic-energy cutoff\s*=\s*([-0-9.]*)'
        outString = AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=energyCutoffRegex)

        return outString

def _grab__volume_atom(oneCalc,ID):
        natom=float(AFLOWpi.aflowlib._grab__natoms(oneCalc,ID))
        vol = AFLOWpi.aflowlib._grab__volume_cell(oneCalc,ID)

        return vol/natom

def _grab__energy_atom(oneCalc,ID):
        natom=float(AFLOWpi.aflowlib._grab__natoms(oneCalc,ID))
        energy = float(AFLOWpi.aflowlib._grab__energy_cell(oneCalc,ID))
        
        return energy/natom

def _grab__natoms(oneCalc,ID):
        outString = ''
        natomsRegex = r'number of atoms[/]*cell\s*=\s*([0-9]*)'
        outString = int(AFLOWpi.aflowlib._outputRegex(oneCalc,ID,natomsRegex))

        return outString

def _grab__nspecies(oneCalc,ID):
        outString = ''
        nspeciesRegex = r'number of atomic types\s*=\s*([0-9]*)'
        outString = AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=nspeciesRegex)

        return outString

def _grab__files(oneCalc,ID):
        outString = ''

        return outString

def _grab__forces(oneCalc,ID):
        return float(AFLOWpi.retr.getForce(oneCalc,ID).split('=')[-1])

def _grab__hubbard_U_values(oneCalc,ID):
        try:
                split_input = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
                system_namelist = split_input['&system']

                u_val_dict={}
                species = AFLOWpi.aflowlib._get_ordering_from_species_list(oneCalc,ID)
                for i in list(system_namelist.keys()):
                        lowered=i.lower()
                        if 'hubbard_u' in lowered:
                                index=int(re.findall('hubbard_u\((\d+)\)',lowered)[-1])
                                u_val_dict[species[index-1]]=float(split_input['&system'][i])

                return_str=''
                for i in species_list:
                        return_str+=str(u_val_dict[i])
                
                return u_val_dict
        except Exception as e:
                pass

def _grab__bravais_lattice_final(oneCalc,ID):
        alat,cellParamMatrix = AFLOWpi.retr._getCellParams(oneCalc,ID)
        cell_vecs=cellParamMatrix*alat
        ibrav = int(AFLOWpi.retr.getIbravFromVectors(cell_vecs))

        return ibrav

def _grab__bravais_lattice_initial(oneCalc,ID):
        cell_vecs=AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        ibrav = int(AFLOWpi.retr.getIbravFromVectors(cell_vecs))

        return ibrav


def _grab__initial_geometry(oneCalc,ID):
        try:
                outString = ''
                ibravRegex = r'bravais-lattice index\s*=\s*([-0-9]*)'
                ibrav = int(AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=ibravRegex)) 
                cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        
                geoString = AFLOWpi.retr.free2abc(cellParamMatrix,cosine=False,)
                commaSplit = geoString.split(',')

                outString = ','.join([item.split('=')[1] for item in commaSplit])

                return outString

        except Exception as e:
                AFLOWpi.run._fancy_error_log(e)
                raise SystemExit

def _grab__final_geometry(oneCalc,ID):
        alat,cellParamMatrix = AFLOWpi.retr._getCellParams(oneCalc,ID)
        cell_vecs=cellParamMatrix*alat

        params = AFLOWpi.retr.free2abc(cell_vecs,cosine=False,string=False)
        out_str=''
        for i in list(params):
                out_str+=i

        return  out_str

def _grab__kpoints(oneCalc,ID):
        outString = ''

        return outString

def _grab__ldaul_TLUJ(oneCalc,ID):
        outString = ''

        return outString


def _grab__initial_lattice_vectors(oneCalc,ID):
        cell_vecs=AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        as_string = AFLOWpi.retr._cellMatrixToString(cell_vecs)
        flattened_str = AFLOWpi.aflowlib._flatten_tensor(as_string,rank=2)

        return flattened_str

def _grab__final_lattice_vectors(oneCalc,ID):
        alat,cellParamMatrix = AFLOWpi.retr._getCellParams(oneCalc,ID)
        cell_vecs=cellParamMatrix*alat
        as_string = AFLOWpi.retr._cellMatrixToString(cell_vecs)
        flattened_str = AFLOWpi.aflowlib._flatten_tensor(as_string,rank=2)

        return flattened_str

def _grab__initial_positions_cartesian(oneCalc,ID):
        pos = AFLOWpi.retr.getPositionsFromInput(oneCalc,ID)
        cell_vecs=AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        pos_cart = AFLOWpi.retr._convertCartesian(pos,cell_vecs)

        as_string = AFLOWpi.retr._cellMatrixToString(pos_cart)
        flattened_str = AFLOWpi.aflowlib._flatten_tensor(as_string,rank=2)

        return flattened_str

def _grab__final_positions_fractional(oneCalc,ID):
        pos = AFLOWpi.retr.getPositionsFromOutput(oneCalc,ID)
        as_string = AFLOWpi.retr._cellMatrixToString(pos)
        flattened_str = AFLOWpi.aflowlib._flatten_tensor(as_string,rank=2)

        return flattened_str

def _grab__initial_positions_fractional(oneCalc,ID):
        pos = AFLOWpi.retr.getPositionsFromInput(oneCalc,ID)
        as_string = AFLOWpi.retr._cellMatrixToString(pos)
        flattened_str = AFLOWpi.aflowlib._flatten_tensor(as_string,rank=2)

        return flattened_str

def _grab__final_positions_cartesian(oneCalc,ID):
        pos = AFLOWpi.retr.getPositionsFromOutput(oneCalc,ID)

        alat,cellParamMatrix = AFLOWpi.retr._getCellParams(oneCalc,ID)
        cell_vecs=cellParamMatrix*alat

        pos_cart = AFLOWpi.retr._convertCartesian(pos,cell_vecs)

        as_string = AFLOWpi.retr._cellMatrixToString(pos_cart)
        flattened_str = AFLOWpi.aflowlib._flatten_tensor(as_string,rank=2)

        return flattened_str




def _grab__pressure(oneCalc,ID):
        outString = ''
        pressureRegex = r'\(kbar\)\s*P=\s*([-0-9.]*)'
        outString = AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=pressureRegex)
        return outString

def _grab__prototype(oneCalc,ID):
        outString = ''

        return outString

def _grab__spacegroup_orig(oneCalc,ID):
        outString = ''

        return outString

def _grab__species(oneCalc,ID):
        atom_list = AFLOWpi.retr._getAtomNum(oneCalc['_AFLOWPI_INPUT_'],strip=True)
        species = list(atom_list.keys())

        return species

def _grab__composition(oneCalc,ID):
        atom_list = AFLOWpi.retr._getAtomNum(oneCalc['_AFLOWPI_INPUT_'],strip=True)
        
        return list(OrderedDict(sorted(list(atom_list.items()), key=lambda key_value: int(key_value[0].split('_')[1]))).keys())

def _grab__species_pp(oneCalc,ID):
        outString = ''

        return outString

def _grab__species_pp_version(oneCalc,ID):
        outString = ''

        return outString

def _grab__spin_cell(oneCalc,ID):
        outString = ''

        return outString

def _grab__spinD(oneCalc,ID):
        outString = ''

        return outString

def _grab__volume_cell(oneCalc,ID):
        outString = ''
        volRegex = r'unit-cell volume\s*=\s*([0-9.]*)'
        outString = float(AFLOWpi.aflowlib._outputRegex(oneCalc,ID,regex=volRegex))

        return outString

def _grab_file_from_prefix(oneCalc,ID,prefix,extension='pdf'):
        try:
                fileplot = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_*%s*.%s'%(prefix,ID,extension))[-1]
                filename = fileplot.split('/')[-1]
                file_directory = fileplot.split('/')[-2]
                file_project = oneCalc['PROJECT']
                file_set = oneCalc['SET']
#               fileplot_rel = os.path.join(file_directory,filename)
                fileplot_rel=filename
                return fileplot_rel

        except Exception as e:
                return None     

def _grab__supercell_dim(oneCalc,ID):
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_fd.in'%ID)) as fd_fo:
                fd_input_str = fd_fo.read()
        nrx1 = int(AFLOWpi.retr._splitInput(fd_input_str)['&inputfd']['nrx1'])
        nrx2 = int(AFLOWpi.retr._splitInput(fd_input_str)['&inputfd']['nrx2'])
        nrx3 = int(AFLOWpi.retr._splitInput(fd_input_str)['&inputfd']['nrx3'])
        out_str='%s,%s,%s' % (nrx1,nrx2,nrx3)


        return out_str

def _grab__file_rham(oneCalc,ID):
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'RHAM',extension='xml')

def _grab__plot_epsilon(oneCalc,ID):
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'EPSILON')
        
def _grab__plot_sigma_seebeck(oneCalc,ID):
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'SIGMA_SEEBECK')

def _grab__plot_conductivity(oneCalc,ID):
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'CONDUCTION')

def _grab__plot_ZetaT(oneCalc,ID): 
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'ZETAT')

def _grab__plot_seebeck(oneCalc,ID):
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'SEEBECK')

def _grab__plot_kappa(oneCalc,ID): 
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'KAPPA')

def _grab__plot_phonon(oneCalc,ID): 
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'PHONON')

def _grab__plot_pdos(oneCalc,ID): 
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'PDOS')

def _grab__plot_dos(oneCalc,ID): 
        return AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'DOS')

def _grab__plot_bands(oneCalc,ID): 
        bands_pdos = AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'BANDPDOS')
        if bands_pdos!=None:
                return bands_pdos
        else:
                bands_dos = AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'BANDDOS')
                if bands_dos!=None:
                        return bands_dos
                else:
                        bands_only = AFLOWpi.aflowlib._grab_file_from_prefix(oneCalc,ID,'BAND')
                        return bands_only

def _grab__epsilon_inf(oneCalc,ID):
        eps_tensor=AFLOWpi.run._pull_eps_out(oneCalc,ID)
        flattened = AFLOWpi.aflowlib._string_to_tensor(eps_tensor,rank=2)
        return  flattened




def _grab__born_eff_charge(oneCalc,ID):
        born_tensor = AFLOWpi.run._pull_born_out(oneCalc,ID)
        atom_list = AFLOWpi.retr._getAtomNum(oneCalc['_AFLOWPI_INPUT_'],strip=False)
        species = list(atom_list.keys())

        by_species=born_tensor.split('\n')
        tensor_list=[]
        by_species_tensor_list=[]
        for i in range(len(by_species)):

                if i%4:
                        tensor_list.append(by_species[i])
                else:
                        if len(tensor_list)==3:
                                by_species_str='\n'.join(tensor_list)
                                by_species_tensor_list.append(AFLOWpi.aflowlib._string_to_tensor(by_species_str,rank=2))
                                
                        tensor_list=[]

        return by_species_tensor_list

def _get_ordering_from_species_list(oneCalc,ID):
        split_input = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
        species_string = split_input['ATOMIC_SPECIES']['__content__']

        species_list = [i.split()[0] for i in species_string.split('\n') if len(i.strip())]

        return species_list 

def _flatten_tensor(string,rank=2):
        if rank==2:
                flattened=';'.join([','.join(i.split()) for i in string.split('\n') if len(i.strip())])

                return flattened

def _string_to_tensor(string,rank=2):
        if rank==2:
                flattened='\n'.join([' '.join(i.split()) for i in string.split('\n') if len(i.strip())])

                return flattened

def _grab__gap_type(oneCalc,ID):
        return AFLOWpi.retr.gap_type(oneCalc,ID)

def _grab__gap_size(oneCalc,ID):
        return AFLOWpi.retr.gap_size(oneCalc,ID)

def _grab__fermi_level(oneCalc,ID):
        return float(AFLOWpi.retr._getEfermi(oneCalc,ID))


def _write_property_list(property_dict,property_list_type):
        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS
        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS
        unit_dict={}
        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS
        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS        #DO THIS


        AFLOWpi_delimiter = '[AFLOWPI] **************************************************************************************************************************'
        prop_list_start='[%s]START'%property_list_type.upper()
        prop_list_stop='[%s]STOP'%property_list_type.upper()
        output_str=''

        output_str+=AFLOWpi_delimiter+'\n'
        output_str+=prop_list_start+'\n'
        for k,v in list(property_dict.items()):
                if k in list(unit_dict.keys()):
                        unit='(%s)'%unit_dict[k]
                else:
                        unit=''
                #do not write to file if the value is null
                try:
                        if len(v.strip())==0:
                                continue
                except:
                        pass
                if v==None or v=='None':
                        continue
                #add property to output string
                output_str+='%s=%s %s\n'%(k,v,unit)

        output_str+=prop_list_stop+'\n'
        output_str+=AFLOWpi_delimiter+'\n'

        return output_str


def _write_properties_for_one(all_properties,dest):
        out_props={}

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('phonon'))
        out_props['phonon']=list(set(props_list))

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('bands'))
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('dos'))
        out_props['electronic']=list(set(props_list))

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('transport'))
        out_props['transport']=list(set(props_list))

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('initial_scf'))
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('final_scf'))
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('only_scf'))
        out_props['structure']=list(set(props_list))

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('PAO-TB'))
        out_props['tb']=list(set(props_list))

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('calculation'))
        out_props['calculation']=list(set(props_list))

        props_list=[]
        props_list.extend(AFLOWpi.aflowlib._get_properties_list('machine'))
        out_props['machine']=list(set(props_list))
        
        for prop_list_name,prop_list in list(out_props.items()):
                property_dict={}
                for i in prop_list:
                        if i in list(all_properties.keys()):
                                property_dict[i]=all_properties[i]

                prop_out_string = AFLOWpi.aflowlib._write_property_list(property_dict,prop_list_name,)
                with open(os.path.join(dest,'aflowpi.%s.out'%prop_list_name),'w') as ofo:
                        ofo.write(prop_out_string)



def _get_properties_list(calc_type):
        #######################################################################################
        #######################################################################################
        ### Master list of possible properties to look for when exporting data from AFLOWpi ###
        #######################################################################################
        #############!!!!! Only edit this if you want to add properties !!!!!##################
        #######################################################################################
        #######################################################################################

        props={}
        props['phonon']      = ['plot_phonon','supercell_dim','epsilon_inf','born_eff_charge']

        props['bands']       = ['plot_bands']

        props['dos']         = ['plot_dos','plot_pdos','gap_type','gap_size','fermi_level']

        props['transport']   = ['plot_kappa','plot_seebeck','plot_sigma_seebeck',
                                'plot_conductivity','plot_epsilon','plot_ZetaT']

        props['final_scf']   = ['volume_cell','pressure','forces','nspecies','natoms','compound',
                                'energy_cell','final_positions_cartesian','final_geometry',
                                'final_positions_fractional','bravais_lattice_final',
                                'hubbard_U_values','final_lattice_vectors',]

        props['only_scf']    = ['nspecies','natoms','bravais_lattice_initial',
                                'initial_positions_cartesian','final_positions_cartesian',
                                'final_geometry','compound','final_positions_fractional',
                                'energy_cell','pressure','forces','bravais_lattice_final',
                                'hubbard_U_values','final_lattice_vectors',
                                'initial_lattice_vectors','position_labels',]
#                               ,'initial_geometry']

        props['initial_scf'] = ['nspecies','natoms',
                                'initial_positions_cartesian','energy_cell','pressure','forces',
                                'initial_positions_fractional','initial_lattice_vectors',
                                'position_labels','bravais_lattice_initial',]
#                               'initial_geometry']
        props['middle_scf']  = []

#       props['scfuj']       = ['uVal','file_kham','file_rham','final_geometry','Bravais_lattice_final','final_atomic_positions',]

        props['PAO-TB']    = ['plot_TB_bands','plot_TB_dos','plot_TB_pdos','file_rham',]

        props['EXCLUDE']   = []
        props['crawl_min']   = []
        props['elastic']   = []
        props['thermal_relax']=[]
        props['thermal']=[]
        #radial_dist=['plot_RDF']

        props['calculation']=['dft_type','auid','energy_cutoff','code',]

        props['machine']=['hostname','kernel_version','architecture']

        return props[calc_type]

def export(project,set_name='',config='',author='',affiliation='',exclude_steps=[],exclude_properties=[],discard_incomplete=True,write=True):
        
#       #load the last step 
#
        calcs=AFLOWpi.aflowlib.load_last(project,set_name,config=config)
        if discard_incomplete==True:
                longest_workflow=0
                for ID,oneCalc in list(calcs.items()):
                        wf_len = len(oneCalc['_AFLOWPI_WORKFLOW_'])
                        if wf_len>longest_workflow:
                                longest_workflow=wf_len
                calcs_copy=copy.deepcopy(calcs)
                for ID,oneCalc in list(calcs_copy.items()):
                        wf_len = len(oneCalc['_AFLOWPI_WORKFLOW_'])
                        if wf_len<longest_workflow:
                                del calcs[ID]


        #dictionary with all the workflows
        export_dict={}

        workdir = AFLOWpi.prep._ConfigSectionMap('prep','work_dir')
        tree_base = os.path.join(workdir,project,set_name)
        export_base = os.path.join(tree_base,'AFLOWpi','__AFLOWPI_EXPORT_TEMP__',project,set_name)

#       ignore_list = shutil.ignore_patterns('_*.*','*.x','*.UPF','.upf','*.py','*.pyc','*.txt','FD_PHONON','GRID_MIN','*.modes','*.sumpdos','CRASH','*.dat','*_reduced_file_*','*.tot','*.efermi','*.in','*.out','*.log','*.gp','AFLOWpi','Nlm_k_file','Sk','Hk_*')

#        include_list = include_patterns('*.pdf', '*.xml')
#       shutil.copytree(tree_base, export_base, symlinks=False, ignore=include_list)

        #generate the empty AFLOWpi folder inside
        if write==True:
            if set_name!='':
                if not os.path.exists(os.path.dirname(os.path.dirname(export_base))):
                        os.mkdir(os.path.dirname(os.path.dirname(export_base)))
                        os.mkdir(os.path.dirname(export_base))
                        os.mkdir(export_base)
#                       os.mkdir(os.path.join(export_base,'AFLOWpi'))
            else:
                if not os.path.exists(os.path.dirname(export_base)):
                        os.mkdir(export_base)
#                       os.mkdir(os.path.join(export_base,'AFLOWpi'))
                
#       export_AFLOWpi_dir = os.path.join(export_base,'AFLOWpi')
#       if not os.path.exists(export_AFLOWpi_dir):
#               os.mkdir(export_AFLOWpi_dir)

        for ID, oneCalc in list(calcs.items()):
                #get workflow list from each calculation in the set
                workflow=oneCalc['_AFLOWPI_WORKFLOW_']
                scf_workflow_list=[]
                scf_list=[]
                index=0
                for i in range(len(workflow)):
                        #split the workflow for each calculation into subsets
                        #every time the charge density changes. If the charge
                        #density changes without properties being calculated with
                        #it then lump the scf steps that occur in succession together
                        if workflow[i] in ['scf','vcrelax','relax','scfuj']:
                                scf_list.append(workflow[i])
                        else:
                                if len(scf_list)!=0:
                                        scf_workflow_list.extend([scf_list])
                                        scf_list=[]

                                scf_workflow_list[-1].append(workflow[i])

                #sort out the initial and final
                #coords for the workflow
                vc_list = ['scf','vcrelax','relax','scfuj']
                ##############################################################################################3
                ##############################################################################################3
                ##############################################################################################3
                ##DONT FORGET TO FIX            ##DONT FORGET TO FIX            ##DONT FORGET TO FIX
                ##DONT FORGET TO FIX            ##DONT FORGET TO FIX            ##DONT FORGET TO FIX
#               scf_workflow_list=scf_workflow_list[-1:]
                ##DONT FORGET TO FXI            ##DONT FORGET TO FIX            ##DONT FORGET TO FIX
                ##DONT FORGET TO FIX            ##DONT FORGET TO FIX            ##DONT FORGET TO FIX
                ##############################################################################################3         
                ##############################################################################################3
                ##############################################################################################3
                for i in range(len(scf_workflow_list)):
                        attribute_list=[]
                        att_dict={}
                        for j in range(len(scf_workflow_list[i])):
                                index+=1
                                #exclude properties or steps if the user had them in the exclude list(s)
                                if scf_workflow_list[i][j] in exclude_properties or index in exclude_steps:
                                        scf_workflow_list[i][j]='EXCLUDE'

                                #set the initial and final scf calcs in the workflows if
                                #there's only one mark it as the only scf type calculation
                                if scf_workflow_list[i][j] in vc_list:
                                        if j==0:
                                                initial_ID = '%s_%02d'%(oneCalc['_AFLOWPI_PREFIX_'][1:-3],index)
                                                scf_workflow_list[i][j]='initial_scf'
                                                if len(scf_workflow_list[i])>1:
                                                        try:
                                                                if scf_workflow_list[i][j+1] not in vc_list:
                                                                        #first and only scf type
                                                                        scf_workflow_list[i][j]='only_scf'
                                                                else:
                                                                        #there are more scf type
                                                                        #so it's not the only one
                                                                        pass
                                                        except:
                                                                pass

                                        if scf_workflow_list[i][j] in vc_list and j!=0:
                                                #if it's a relax and it's not the first entry check
                                                #to see if the next entry is also some kind of scf
                                                try:
                                                        if scf_workflow_list[i][j+1] not in vc_list:
                                                                #if this is the last scf type in the list the
                                                                #list then we do final cell parameters..etc
                                                                scf_workflow_list[i][j]='final_scf'
                                                        else:
                                                                scf_workflow_list[i][j]='middle_scf'
                                                except:
                                                        #must be the only entry in this list 
                                                        #so it's both initial and final
                                                        scf_workflow_list[i][j]='only_scf'

                                #gather an attribute list for the type of 
                                #calculation done step done on this step.
                                attribute_list=AFLOWpi.aflowlib._get_properties_list(scf_workflow_list[i][j])
                                index_ID = '%s_%02d'%(oneCalc['_AFLOWPI_PREFIX_'][1:-3],index)
                                new_hash = oneCalc['_AFLOWPI_PREFIX_'][1:-3]

                                for att in attribute_list:
                                        val = AFLOWpi.aflowlib._get_attribute(oneCalc,index_ID,att)
                                        att_dict[att]=val

                        #get new hash by combining input of first calc in the subset
                        #with the output of it. 
                        hash_string =  AFLOWpi.retr._getInputFileString(oneCalc,initial_ID)
                        hash_string += AFLOWpi.retr._getOutputString(oneCalc,initial_ID)
#                       new_hash = AFLOWpi.prep._hash64String(hash_string)
                        export_dict[new_hash] = att_dict



        # for ID,oneCalc in export_dict.iteritems():
        #       print 
        #       print
        #       print ID
        #       for k,v in oneCalc.iteritems():
        #               print k,v

        if write==False:
            return export_dict

        folder_index=0
        #remake the 
        for k,v in list(export_dict.items()):
                folder_index+=1
                one_dir_str = "%s_%s_%04d"%(project,set_name,folder_index)
                one_dir = os.path.join(export_base,one_dir_str)

                try:
                        os.mkdir(one_dir)
                except:
                        pass
                AFLOWpi.aflowlib._write_properties_for_one(v,one_dir)           

                data_dir= os.path.join(tree_base,one_dir_str)
                AFLOWpi.aflowlib._copy_data_dir(data_dir,one_dir)

        
        #dump the dictionary of dictionaries containing 
        #information about the calc set being exported
#       export_dict_filename = os.path.join(export_base,'AFLOWpi','set_data.pkl')
#       with open(export_dict_filename,'wb') as export_dict_file_obj:                   
#               cPickle.dump(export_dict,export_dict_file_obj)    



        #create the export tar file and dump the contents
        #of the temp export directory in it.
        if set_name=='':
                tar_name='%s'%(project)
        else:
                tar_name='%s_%s'%(project,set_name)

        prev_dir=os.chdir
        os.chdir(os.path.join(tree_base,'AFLOWpi'))

        #archive the export in a tar for transport and storage
        shutil.make_archive(tar_name,'gztar','./__AFLOWPI_EXPORT_TEMP__/%s'%project)

        #remove the temp export directory
        temp_base = os.path.join(tree_base,'AFLOWpi','__AFLOWPI_EXPORT_TEMP__')
        shutil.rmtree(temp_base)

        try:
                os.chdir(prev_dir)
        except:
                pass

                
def _get_attribute(oneCalc,ID,attribute):
        try:
                val=eval('AFLOWpi.aflowlib._grab__%s(oneCalc,ID)'%attribute)
                return val
        except Exception as e:
                return None

def load_last(project,set_name='',config=''):

        calc_list=AFLOWpi.retr.checkStatus(project,SET=set_name,config=config,step=0,negate_status=True)
        last_dict={}
        for i in range(len(calc_list)):
                for k,v in list(calc_list[i].items()):
                        last_dict[v["_AFLOWPI_PREFIX_"][1:].split("_")[0]]=v

        return last_dict

def _copy_data_dir(data_dir,dest):
        copy_list=[]
        copy_list.extend(glob.glob(os.path.abspath(data_dir)+'/*.pdf'))
        copy_list.extend(glob.glob(os.path.abspath(data_dir)+'/*.in'))
        copy_list.extend(glob.glob(os.path.abspath(data_dir)+'/*.out'))
#       copy_list.extend(glob.glob(os.path.abspath(data_dir)+'/*.xml'))

        for file_name in copy_list:
                try:            
                        shutil.copy(file_name,dest)
                except Exception as e:
                        print(e)
                        pass

# def include_patterns(*patterns):
#     """Factory function that can be used with copytree() ignore parameter.

#     Arguments define a sequence of glob-style patterns
#     that are used to specify what files to NOT ignore.
#     Creates and returns a function that determines this for each directory
#     in the file hierarchy rooted at the source directory when used with
#     shutil.copytree().
#     """
#     def _ignore_patterns(path, names):
#         keep = set(name for pattern in patterns
#                             for name in fnmatch.filter(names, pattern))
#         ignore = set(name for name in names
#                         if name not in keep and not os.path.isdir(os.path.join(path, name)))
#         return ignore
#     return _ignore_patterns
