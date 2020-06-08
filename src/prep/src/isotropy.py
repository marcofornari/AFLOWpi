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
import numpy as np
import os
import subprocess
import io
import re
import copy
import scipy


np.set_printoptions(precision=4, threshold=200, edgeitems=200, linewidth=250, suppress=True)                   
class isotropy():

    def __init__(self):


        self.iso_input = ''
        self.qe_basis  = ''
        self.pos_labels= ''
        self.orig_pos  = ''
        self.iso_basis = ''
        self.qe_pos    = ''

        self.conv_a    = ''
        self.conv_c    = ''
        self.conv_b    = ''
        self.conv_alpha= ''
        self.conv_beta = ''
        self.conv_gamma= ''
        self.origin    = ''
        self.cif       = ''

    def qe_input(self,input_str,accuracy=0.001):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()

        
        self.input     =  AFLOWpi.prep._removeComments(input_str)
        self.accuracy  = accuracy

        self.output    = self.__get_isotropy_output(qe_output=False)
        self.cif       = self.qe2cif()
        self.sgn       = self.__get_sg_num()
        self.ibrav     = self.__ibrav_from_sg_number()
        self.iso_pr_car= ''
        
#        self.iso_basis = self.__get_iso_basis()
#        self.iso_pr_car= self.__get_iso_cart()

        input_dict = AFLOWpi.retr._splitInput(self.input)

    def qe_output(self,input_str,accuracy=0.001):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()


        self.input     =  AFLOWpi.prep._removeComments(input_str)
        self.accuracy  = accuracy

        self.output    = self.__get_isotropy_output(qe_output=True)
        self.cif       = self.qe2cif()
        self.sgn       = self.__get_sg_num()
        
    
        self.ibrav     = self.__ibrav_from_sg_number()
        self.iso_pr_car= ''

        self.iso_basis = self.__get_iso_basis()
        self.iso_pr_car= self.__get_iso_cart()

    def cif_input(self,input_str):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()

        self.accuracy  = 0.0001
        self.cif       = input_str
        self.output    = input_str
        self.input     = self.__skeleton_input()
        self.accuracy  = 0.0001
        self.sgn       = self.__get_sg_num()
        
        self.ibrav     = self.__ibrav_from_sg_number()
        self.origin    = self.__conv_from_cif()        

        self.iso_pr_car= ''


    def __skeleton_input(self):
        skel_in = '''
&control
/
&system
ibrav=-1
nat=0
ntyp=0
ecutwfc=40.0
ecutrho=400.0
/
&electrons
/
&ions
/
&cell
/
ATOMIC_SPECIES

K_POINTS {automatic}
'''
        return skel_in

    def qe2cif(self):
        '''
        Converts Quantum Espresso input into CIF
        
        Returns: (str) CIF file generated from Quantum Espresso input
        '''
        # we need to take the wyckoff positions and match
        # the ATOMIC_SPECIES labels in the qe input
        # wyc_full = re.findall('Atomic positions in terms of a,b,c:\n((?:.*\n)+)------',
        #                       self.output)[0]
        # each_wyc = re.findall('Wyckoff position.*\n(?:(\s*\d+\s*[\d.\s]+\n)+)',wyc_full)

        # wyc_first_lab_ind=[]
        # for i in each_wyc:
        #     wyc_first_lab_ind.append([int(x.split()[0])-1 for x in i.split('\n') if len(x.strip())!=0][0])

        # cif_file = re.findall('# CIF file.*\n(?:.*\n)*',self.output)[0]
        # re_atom_pos=re.compile(r'(_atom_site_label.*\n(?:(?:[A-Za-z_\s])*\n))((?:.*\n)*)')

        # atom_pos = re_atom_pos.findall(cif_file)[0]
        # print(atom_pos)
        # raise System
        # '''positions from cif'''
        # loop_list = [x.strip() for x in  atom_pos[0].split('\n') if len(x.strip())!=0]
        # spec_lab_index =  loop_list.index('_atom_site_type_symbol')

        # positions = [list(map(str.strip,x.split())) for x in  atom_pos[1].split('\n') if (len(x.strip())!=0 and len(x.split())==len(loop_list))]

        # mod_pos=[]
        # for i in range(len(positions)):
        #     positions[i][spec_lab_index]=self.pos_labels[wyc_first_lab_ind[i]]
        #     mod_pos.append(' '.join(['%8.8s'%x for x in positions[i]]))


        # atom_pos[0].split('\n')
        # cif_file_temp_split = [x for x in cif_file.split('\n') if len(x.strip())!=0]
        # cif_file_temp_split[-len(positions):] = mod_pos

        # cif_file = '\n'.join(cif_file_temp_split)+'\n'





        return self.output

    def __generate_isotropy_input_from_qe_data(self,output=False):
        if output:

            cm_string = AFLOWpi.qe.regex.cell_parameters(self.input,'content') 

            alatSearch = re.compile(r'(?:CELL_PARAMETERS)\s*.*alat[\D]*([0-9.]*)',re.M)
            try:
                alat = float(alatSearch.findall(self.input)[-1])
            except:
                alat=1.0

            cell_matrix = AFLOWpi.retr._cellStringToMatrix(cm_string)

            cell_matrix*= alat*0.529177249
            cm_string = AFLOWpi.retr._cellMatrixToString(cell_matrix)

            pos_with_labs = AFLOWpi.qe.regex.atomic_positions(self.input,'content') 

            labels=[]
            positions = []
            
            for i in pos_with_labs.split('\n'):
                split_pos = i.split()
                try:
                    labels.append(split_pos[0])
                    split_pos[1]=str(float(split_pos[1])%1.0)
                    split_pos[2]=str(float(split_pos[2])%1.0)
                    split_pos[3]=str(float(split_pos[3])%1.0)
                    positions.append(' '.join(split_pos[1:4]))
                except:
                    pass
            positions='\n'.join(positions)

        else:
            cell_matrix = AFLOWpi.retr.getCellMatrixFromInput(self.input,string=False)*0.529177249
            positions = AFLOWpi.retr._getPositions(self.input,matrix=True)

            try:
                mod = AFLOWpi.qe.regex.atomic_positions(self.input,'modifier')
            except: mod=""

            if mod.lower() in ["angstrom","bohr"]:
                if mod.lower()=="bohr":
                    positions*=0.529177249                
                positions = np.linalg.inv(cell_matrix).T.dot(positions.T).T

            positions = AFLOWpi.retr._cellMatrixToString(positions,indent=False)
        
            labels = AFLOWpi.retr._getPosLabels(self.input)
            self.orig_pos = AFLOWpi.retr._getPositions(self.input,matrix=True)
        
            cm_string = AFLOWpi.retr._cellMatrixToString(cell_matrix,indent=False)
            a,b,c,alpha,beta,gamma = AFLOWpi.retr.free2abc(cell_matrix,cosine=False,bohr=False,string=False)

        self.pos_labels=labels
        num_atoms=len(labels)
        spec = list(set(labels))
        
        label_index=dict([[i[1],i[0]+1] for i in enumerate(spec)])

        in_list={}

        isotropy_input_str='input file for isotropy generated from pwscf input by AFLOWpi\n'    
        isotropy_input_str+='%s\n'%self.accuracy
        isotropy_input_str+='1\n'
        isotropy_input_str+=cm_string

        isotropy_input_str+='2\n'
        isotropy_input_str+='P'+'\n'
        isotropy_input_str+='%s\n'%num_atoms
        isotropy_input_str+=' '.join([str(j) for j in labels])+'\n'

        isotropy_input_str+=positions

        self.iso_input = isotropy_input_str

        return isotropy_input_str

    def __get_isotropy_output(self,qe_output=False):


        # centering='P'

        # in_str = self.__generate_isotropy_input_from_qe_data(output=qe_output)
        # ISODATA = os.path.join(AFLOWpi.__path__[0],'ISOTROPY')
        # os.putenv('ISODATA',ISODATA+'/') 
        # findsym_path = os.path.join(ISODATA,'findsym')


        # try:
        #     find_sym_process = subprocess.Popen(findsym_path,stdin=subprocess.PIPE,stdout=subprocess.PIPE,)
        #     output = find_sym_process.communicate(input=in_str)[0]
        #     self.output=output
        # except:
        #     print((find_sym_process.returncode))

        return  AFLOWpi.retr._get_cif_aflow(self.input,thresh=self.accuracy)

    def __get_sg_num(self):
        try:
            sg_info = re.findall('Space Group\s*([0-9]*)\s*([\w-]*)\s*([\w-]*)',self.output)[0]
            self.sgn = int(sg_info[0])
            return self.sgn
        except:
            pass
        try:
                sg_info = re.findall('_symmetry_Int_Tables_number\s*(\d+)',self.cif)[0]
                self.sgn = int(sg_info)
                return self.sgn
        except:
            pass

        try:
            sg_nam = re.compile("_symmetry_space_group_name_H-M\s*(.*)\n")
            name=sg_nam.findall(self.cif)[0].replace("'","").replace("\r","").replace(" ","")

            self.sgn = self.__HM2Num(name)
            return self.sgn

        except Exception as e:
            return 1


    def __get_iso_cart(self):
            search = 'Lattice vectors in cartesian coordinates:\s*\n(.*\n.*\n.*)\n'
            std_prim_basis_str = re.findall(search,self.output)[0]
            self.iso_pr_car    = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)

            return self.iso_pr_car
            
    def __conv_from_cif(self):
        self.conv_a = float(re.findall('_cell_length_a\s*([0-9-.]*)',self.output)[0])
        self.conv_b = float(re.findall('_cell_length_b\s*([0-9-.]*)',self.output)[0])
        self.conv_c = float(re.findall('_cell_length_c\s*([0-9-.]*)',self.output)[0])
        self.conv_alpha = float(re.findall('_cell_angle_alpha\s*([0-9-.]*)',self.output)[0])
        self.conv_beta  = float(re.findall('_cell_angle_beta\s*([0-9-.]*)',self.output)[0])
        self.conv_gamma = float(re.findall('_cell_angle_gamma\s*([0-9-.]*)',self.output)[0])
        try:
            origin = re.findall('Origin at\s*([0-9-.]*)\s*([0-9-.]*)\s*([0-9-.]*)',self.output)[0]
        except:
            origin=[0.0,0.0,0.0]
        return origin

    def __get_iso_basis(self):
        return np.eye(3)

        # self.iso_basis

        # modifier = AFLOWpi.qe.regex.cell_parameters(self.input,return_which='modifier',).lower()

        # origin = self.__conv_from_cif()

        # origin = np.asarray([float(i) for i in origin])
        # self.origin=origin

        # std_prim_basis_str = re.findall('Vectors a,b,c:\s*\n(.*\n.*\n.*)\n',self.output)[0]
        # self.iso_conv = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)

        # input_dict = AFLOWpi.retr._splitInput(self.input)
        # if 'CELL_PARAMETERS' not in input_dict:
        #     prim_in = AFLOWpi.retr.getCellMatrixFromInput(self.input)
        #     try:
        #         prim_in=AFLOWpi.retr._cellStringToMatrix(prim_in)
        #     except:
        #         pass
        # else:
        #     if input_dict['CELL_PARAMETERS']["__content__"]=='':
        #         try:
        #             prim_in = AFLOWpi.retr.getCellMatrixFromInput(self.input)
        #         except:
        #             pass
        #     try:
        #         prim_in=AFLOWpi.retr._cellStringToMatrix(input_dict['CELL_PARAMETERS']['__content__'])
        #     except Exception as e:
        #         print(e)
        #         raise SystemExit

        # self.iso_basis=prim_in

        # a=np.array([[self.conv_a,],
        #                [self.conv_b,],
        #                [self.conv_c,],])

        # a=np.abs((self.iso_conv).dot(self.iso_basis))

        # prim_abc = re.findall('Lattice parameters, a,b,c,alpha,beta,gamma.*\n(.*)\n',self.output)[0].split()
        # prim_abc=list(map(float,prim_abc))

        # return 


    def __ibrav_from_sg_number(self):
        if self.sgn in [195,198,200,201,205,207,208,212,213,215,218,221,222,223,224]:
            return 1
        elif self.sgn in [196,202,203,209,210,216,219,225,226,227,228]:
            return 2
        elif self.sgn in [197,199,204,206,211,214,217,220,229,230,]:
            return 3
        elif self.sgn in [168,169,170,171,172,173,174,175,176,177,178,179,
                     180,181,182,183,184,185,186,187,188,189,190,191,
                     192,193,194]:
            return 4


        elif self.sgn in [143,144,145,147,149,150,151,152,153,154,
                     156,157,158,159,162,163,164,165,]:
            return 4

        elif self.sgn in [146,148,155,160,161,166,167]:
            return 5
                     
        elif self.sgn in [75,76,77,78,81,83,84,85,86,89,91,92,93,94,
                     95,96,99,100,101,102,103,104,105,106,111,112,
                     113,114,115,116,117,118,123,124,125,126,127,
                     128,129,130,131,132,133,134,135,136,137,138]:
            return 6
        elif self.sgn in [79,80,82,87,88,90,97,98,107,108,109,110,119,
                     120,121,122,139,140,141,142]:
            return 7

        elif self.sgn in [16,17,18,19,25,26,27,28,29,30,31,32,33,34,47,
                     48,49,50,51,52,53,54,55,56,57,58,59,60,61,62]:
            return 8
        elif self.sgn in [38,39,40,41,20,21,35,36,37,63,64,65,66,67,68]:
            return 9

        elif self.sgn in [22,42,43,69,70]:
            return 10
        elif self.sgn in [23,24,44,45,46,71,72,73,74]:
            return 11
        elif self.sgn in [3,4,6,7,10,11,13,14]:
            return 12
        elif self.sgn in [5,8,9,12,15]:
            return 13
        elif self.sgn in [1,2]:
            return 14

        else:
            print(('SG num not found:',self.sgn))
            raise SystemExit



    def convert(self,ibrav=True,thresh=0.01):
        inputString=''

        for thresh in [0.1,0.01,0.1]:
            inputString_new= self.cif2qe(thresh=thresh)
            isotropy_chk    = AFLOWpi.prep.isotropy()
            isotropy_chk.qe_input(inputString_new,accuracy=self.accuracy)


            inputString=inputString_new
        
            if self.sgn!=isotropy_chk.sgn:
                print(("warning "*10))
                print(("warning "*10))
                print(("warning "*10))
                print(inputString)
                print(("warning "*10))
                print(("warning "*10))
                print(("warning "*10))
            break


        return inputString


    def cif2qe(self,thresh=0.01):

        input_dict = AFLOWpi.retr._splitInput(self.input)

        ins= self.cif.lower()
        '''grab the symmetry operations from the text in the cif'''
        re_symops=re.compile(r'''_space_group_symop_operation_xyz\s*\n((?:\s*\d*\s*[\'\"xyz\s,\d\/+-]+\n)+)\s*''',re.M)

        symops = [x for x in re_symops.findall(ins)[0].split('\n') if (len(x.strip())!=0 )]

        re_sym_ops_remove_numbering = re.compile('\d*\s+(.*)')
        for i in  range(len(symops)):
            symops[i] = re_sym_ops_remove_numbering.findall(symops[i])[0]
            symops[i] = symops[i].replace("'","")
            symops[i] = symops[i].replace('"','')


        symops_aux=[]
        for i in  range(len(symops)):
            sym_aux_temp = [x.strip().lower() for x in symops[i].split(',') ]
            if len(sym_aux_temp)!=0:
                symops_aux.append(sym_aux_temp)

        '''prepare to get the symOps'''
        symops=symops_aux
        ident = np.identity(3,dtype=np.float64)

        shift,operations = AFLOWpi.prep._process_cif_sym(symops)
        re_atom_pos=re.compile(r'(_atom_site_label.*\n(?:(?:[a-z_\s])*\n))((?:.*\n)*)')
        atom_pos=re_atom_pos.findall(ins)[0]
        '''positions from cif'''
        loop_list = [x.strip() for x in  atom_pos[0].split('\n') if len(x.strip())!=0]
        
        spec_lab =  loop_list.index('_atom_site_type_symbol')
        x_loc =  loop_list.index('_atom_site_fract_x')
        y_loc =  loop_list.index('_atom_site_fract_y')
        z_loc =  loop_list.index('_atom_site_fract_z')
        try:
            occup =  loop_list.index('_atom_site_occupancy')
        except:
            pass
        '''get positions from input'''
        positions = [list(map(str.strip,x.split())) for x in  atom_pos[1].split('\n') if (len(x.strip())!=0 and len(x.split())==len(loop_list))]

        pos_array=np.zeros((len(positions),3))
        
        labels=[]
        for i in range(len(positions)):
            pos_array[i][0] = float(positions[i][x_loc])
            pos_array[i][1] = float(positions[i][y_loc])
            pos_array[i][2] = float(positions[i][z_loc])
            try:
                occupancy = float(positions[i][occup])
                if occupancy!=1.0:
                    raise IndexError
            except IndexError:
                return
            except:
                pass

            element_list=["H", "He","Li","Be","B", "C", "N", "O", "F", "Ne","Na","Mg","Al",
                          "Si","P", "S", "Cl","Ar","K", "Ca","Sc","Ti","V", "Cr","Mn","Fe",
                          "As","As","Se","Br","Kr","Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru",
                          "Sb","Te","I", "Xe","Cs","Ba","La","Ce","Pr","Nd","Sm","Eu","Gd",
                          "Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi",
                          "Co","Ni","Cu","Zn","Ga","Ge","Rh","Pd","Ag","Cd","In","Sn","Tb",
                          "Dy","Ho","Er","Tm","Yb","Po","At","Th","Pa","U","Am",]


            pos_lab_strip = positions[i][spec_lab].strip('0123456789').title()

            

            try:
                if pos_lab_strip not in element_list:
                    raise IndexError


            except IndexError:                
                return

            
            temp_pos_lab = positions[i][spec_lab].title()

            labels.append(temp_pos_lab)        



        '''shift positions to between 0.0 and 1.0'''
        pos_array%=1.0

        '''if positions is 1/3 or 2/3 use more precision on position'''
        pos_array[np.where(np.isclose(pos_array-(1.0/3.0),0.0,rtol=1e-04, atol=1e-05))]=1.0/3.0
        pos_array[np.where(np.isclose(pos_array-(2.0/3.0),0.0,rtol=1e-04, atol=1e-05))]=2.0/3.0

        '''there will by natm*numSymOps'''
        labels = labels*len(operations)
        all_eq_pos=np.zeros((pos_array.shape[0]*operations.shape[0],3))

        '''duplicate atoms by symmetry operations'''
        for i in range(operations.shape[0]):

            temp_pos=np.copy(pos_array)
            temp_pos   = temp_pos.dot(operations[i])
            '''translation symmetry'''
            temp_pos[:,0]+=shift[i][0]
            temp_pos[:,1]+=shift[i][1]
            temp_pos[:,2]+=shift[i][2]

            '''add atoms generated by symmetry operations to the list of atoms'''
            all_eq_pos[i*pos_array.shape[0]:(i+1)*pos_array.shape[0]]=temp_pos



        '''grab the conventional -> primitive conversion matrix for this ibrav'''
        if self.ibrav in [4,5]:
            t_ibrav=1
        else:
            t_ibrav=self.ibrav

        '''conventional to qe convention primitive lattice vec'''
        convert = AFLOWpi.retr.abc2free(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,
                                        gamma=90.0,ibrav=t_ibrav,returnString=False)

        '''#####################'''
        '''### SPECIAL CASES ###'''
        '''#####################'''

        '''trigonal P hex to rho vecs''' 
        if self.ibrav==5:

            '''if defined with hexagonal axes..transform to rhom'''
            if self.conv_gamma==120.0 or self.conv_beta==120.0 or self.conv_alpha==120.0:
                beta=np.sqrt(3.0+(self.conv_c/self.conv_a)**2.0)
                self.conv_a = self.conv_a*beta/3.0
                self.conv_b = self.conv_a
                self.conv_c = self.conv_a
                self.conv_alpha =  2.0*np.arcsin((3.0/(2.0*beta)))*(180.0/np.pi)
                self.conv_beta  = self.conv_alpha
                self.conv_gamma = self.conv_alpha

                '''hex to rho'''
                convert = np.array([[-1.,  1., -0.,],
                                    [ 1.,  0., -1.,],
                                    [ 1.,  1.,  1.,],])

                convert=np.linalg.inv(np.around(convert,decimals=1))

        """A base centered to C"""
        if self.sgn in [38,39,40,41]:
            all_eq_pos=all_eq_pos[:,[2,1,0]]
            all_eq_pos[:,0]*=-1.0
            temp_c = self.conv_c
            self.conv_c=self.conv_a
            self.conv_a=temp_c

            convert = AFLOWpi.retr.abc2free(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,
                                            gamma=self.conv_gamma,ibrav=self.ibrav,returnString=False)



        '''make A<B<C'''
        if self.ibrav in [14]:           
            cell_vec=AFLOWpi.retr.abc2free(self.conv_a,self.conv_b,self.conv_c,self.conv_alpha,
                                           self.conv_beta,self.conv_gamma,self.ibrav,returnString=False)

            prim_a = (np.sqrt(cell_vec[0].dot(cell_vec[0].T))).getA()[0][0]
            prim_b = (np.sqrt(cell_vec[1].dot(cell_vec[1].T))).getA()[0][0]
            prim_c = (np.sqrt(cell_vec[2].dot(cell_vec[2].T))).getA()[0][0]

            sizes = np.array([prim_a,prim_b,prim_c])
            order = np.argsort(sizes)

            new_cell_vec =  cell_vec[order]

            prim_a = (np.sqrt(new_cell_vec[0].dot(new_cell_vec[0].T))).getA()[0][0]
            prim_b = (np.sqrt(new_cell_vec[1].dot(new_cell_vec[1].T))).getA()[0][0]
            prim_c = (np.sqrt(new_cell_vec[2].dot(new_cell_vec[2].T))).getA()[0][0]

            self.conv_alpha = np.arccos(new_cell_vec[1].dot(new_cell_vec[2].T).getA()[0][0]/(prim_b*prim_c))*180/np.pi
            self.conv_beta  = np.arccos(new_cell_vec[1].dot(new_cell_vec[2].T).getA()[0][0]/(prim_a*prim_c))*180/np.pi
            self.conv_gamma = np.arccos(new_cell_vec[0].dot(new_cell_vec[1].T).getA()[0][0]/(prim_a*prim_b))*180/np.pi
            self.conv_a = prim_a
            self.conv_b = prim_b
            self.conv_c = prim_c

            all_eq_pos=all_eq_pos[:,order]


        '''make unique a or b monoclinic into unique c'''
        if self.ibrav in [12,13]:
            if np.isclose(self.conv_beta,self.conv_gamma) and not np.isclose(self.conv_beta,self.conv_alpha):
                all_eq_pos=all_eq_pos[:,[1,2,0]]
                temp_a=self.conv_b
                temp_b=self.conv_c
                temp_c=self.conv_a
                temp_alpha=self.conv_beta
                temp_beta =self.conv_gamma
                temp_gamma=self.conv_alpha
                self.conv_a=temp_a
                self.conv_b=temp_b
                self.conv_c=temp_c
                self.conv_alpha = temp_alpha
                self.conv_beta  = temp_beta
                self.conv_gamma = temp_gamma

            elif np.isclose(self.conv_alpha,self.conv_gamma) and not np.isclose(self.conv_beta,self.conv_gamma):
                all_eq_pos=all_eq_pos[:,[2,0,1]]
                
                temp_a=self.conv_c
                temp_b=self.conv_a
                temp_c=self.conv_b
                temp_alpha=self.conv_gamma
                temp_beta =self.conv_alpha
                temp_gamma=self.conv_beta
                self.conv_a=temp_a
                self.conv_b=temp_b
                self.conv_c=temp_c
                self.conv_alpha = temp_alpha
                self.conv_beta  = temp_beta
                self.conv_gamma = temp_gamma


        if self.ibrav==13:
            self.ibrav=12
            
        if self.ibrav==12:
            convert=np.matrix(((1.0, 0.0,0.0,),
                               (0.0, 1.0, 0.0,),
                               (0.0, 0.0, 1.0,),))
        '''#####################'''
        '''### SPECIAL CASES ###'''
        '''#####################'''


        '''shift origin to zero before transforming'''
        all_eq_pos-=self.origin


        '''remove duplicate atoms'''

        '''shift all atoms back into cell if need be'''
        all_eq_pos%=1.0


        '''order the the axes a<b<c for ortho'''
        if self.ibrav in [8,10,11]:
             sides=np.array([self.conv_a,self.conv_b,self.conv_c])
             order = np.argsort(sides)
             ordered_sides = sides[order]
             self.conv_a = ordered_sides[0]
             self.conv_b = ordered_sides[1]
             self.conv_c = ordered_sides[2]

        elif self.ibrav==9:
            if self.conv_b<self.conv_a:
                temp_a=self.conv_a
                self.conv_a=self.conv_b
                self.conv_b=temp_a            
                order=np.array([1,0,2])
            else:
                order=np.array([0,1,2])

        else:
            order=np.array([0,1,2])

        all_eq_pos=all_eq_pos[:,order]


        '''convert from conventional to primitive lattice coords'''
        try:
            all_eq_pos=(np.linalg.inv(convert.getA()).T.dot(all_eq_pos.T)).T
        except:
            all_eq_pos=(np.linalg.inv(convert).T.dot(all_eq_pos.T)).T

        all_eq_pos%=1.0

        '''to convert from primitive to conventional for distance calc'''
        to_conv = AFLOWpi.retr.abc2free(a=self.conv_a,b=self.conv_b,c=self.conv_c,
                                        alpha=self.conv_alpha,beta=self.conv_beta,
                                        gamma=self.conv_gamma,ibrav=self.ibrav,
                                        returnString=False)


        if self.ibrav==13:
            two_one = np.cos(self.conv_gamma/180.0*np.pi)*self.conv_b
            two_two = np.sin(self.conv_gamma/180.0*np.pi)*self.conv_b

            to_conv=np.matrix(((1.0*self.conv_a, 0.0,0.0,),
                               (0.0, 1.0*self.conv_b, 0.0,),
                               (0.0, 0.0, 1.0*self.conv_c,),))

        if self.ibrav==12:
            to_conv=np.matrix(((1.0*self.conv_a, 0.0,0.0,),
                               (0.0, 1.0*self.conv_b, 0.0,),
                               (0.0, 0.0, 1.0*self.conv_c,),))



        '''reduce equivilent atoms'''
        all_eq_pos,labels_arr = reduce_atoms(all_eq_pos,labels,cell=to_conv,thresh=thresh)


        '''form atomic positions card for QE input'''
        atm_pos_str=""
        for i in range(all_eq_pos.shape[0]):
            atm_pos_str+= ('%3.3s'% labels_arr[i]) +(' % 9.9f % 9.9f % 9.9f '%tuple(all_eq_pos[i].tolist()))+"\n"


        '''assign values of cell to A,B,C,cosAB,cosAC,cosBC in QE input file'''
        input_dict['&system']['ibrav']=self.ibrav

        input_dict['&system']['A']=self.conv_a
        input_dict['&system']['B']=self.conv_b
        input_dict['&system']['C']=self.conv_c

        input_dict['&system']['cosAB']=np.cos(self.conv_gamma/180.0*np.pi)
        input_dict['&system']['cosAC']=np.cos(self.conv_beta/180.0*np.pi)
        input_dict['&system']['cosBC']=np.cos(self.conv_alpha/180.0*np.pi)


        """get rid of cell parameters card if present"""
        try:
            del input_dict['CELL_PARAMETERS']
        except: pass
        
        """assign new number of atoms"""
        input_dict['&system']['nat']=all_eq_pos.shape[0]
        
        '''do ntyp'''
        input_dict['&system']['ntyp']=len(set(labels))

        '''add newly transformed atomic positions to QE convention input'''
        try:
            input_dict['ATOMIC_POSITIONS']['__content__']=atm_pos_str        
        except:
            input_dict['ATOMIC_POSITIONS']={}
            input_dict['ATOMIC_POSITIONS']["__modifier__"]="{crystal}"
            input_dict['ATOMIC_POSITIONS']['__content__']=atm_pos_str        


        # """reorder MP kpoint grid if axes change ordering"""
        # mod = input_dict['K_POINTS']['__modifier__'].strip("{").strip("}").strip("(").strip(")").lower()
        # if mod=="automatic":
        #     kp = np.array(map(int,input_dict['K_POINTS']['__content__'].split()))
        #     cell_vec=AFLOWpi.retr.abc2free(self.conv_a,self.conv_b,self.conv_c,self.conv_alpha,
        #                                    self.conv_beta,self.conv_gamma,self.ibrav,returnString=False)

        #     new_grid = AFLOWpi.prep.getMPGrid(cell_vec/0.529177249,offset=True,string=True)
        #     '''recip lattice'''
        #     cell_vec=np.linalg.inv(cell_vec)
        #     '''prim recip lattice vectors lengths'''

        #     prim_a = (np.sqrt(cell_vec[0].dot(cell_vec[0].T))).getA()[0][0]
        #     prim_b = (np.sqrt(cell_vec[1].dot(cell_vec[1].T))).getA()[0][0]
        #     prim_c = (np.sqrt(cell_vec[2].dot(cell_vec[2].T))).getA()[0][0]

        #     '''find the order from smallest to largest of the recip lattice vec'''
        #     ordering = np.argsort(np.array([prim_a,prim_b,prim_c]))

        #     '''order new kp grid with largest for smallest real space vector'''
        #     kp[:3][ordering] = np.sort(kp[:3])
        #     input_dict['K_POINTS']['__content__'] =  ' '.join(map(str,kp.tolist()))
        #     input_dict['K_POINTS']['__content__'] =  new_grid

        qe_convention_input=AFLOWpi.retr._joinInput(input_dict)
        '''convert to celldm'''
        qe_convention_input = AFLOWpi.prep._transformInput(qe_convention_input)


        return qe_convention_input


    def __HM2Num(self,HM_str):
        HMdict={"P1" : 1,
                "A1" : 1,
                "B1" : 1,
                "C1" : 1,
                "I1" : 1,
                "F1" : 1,
                "P-1" : 2,
                "A-1" : 2,
                "B-1" : 2,
                "C-1" : 2,
                "I-1" : 2,
                "F-1" : 2,
                "P2" : 3,
                "P2" : 3,
                "P2" : 3,
                "P21" : 4,
                "P21" : 4,
                "P21" : 4,
                "C2" : 5,
                "A2" : 5,
                "I2" : 5,
                "A2" : 5,
                "B2" : 5,
                "I2" : 5,
                "B2" : 5,
                "C2" : 5,
                "I2" : 5,
                "Pm" : 6,
                "Pm" : 6,
                "Pm" : 6,
                "Pc" : 7,
                "Pn" : 7,
                "Pa" : 7,
                "Pa" : 7,
                "Pn" : 7,
                "Pb" : 7,
                "Pb" : 7,
                "Pn" : 7,
                "Pc" : 7,
                "Cm" : 8,
                "Am" : 8,
                "Im" : 8,
                "Am" : 8,
                "Bm" : 8,
                "Im" : 8,
                "Bm" : 8,
                "Cm" : 8,
                "Im" : 8,
                "Cc" : 9,
                "An" : 9,
                "Ia" : 9,
                "Aa" : 9,
                "Bn" : 9,
                "Ib" : 9,
                "Bb" : 9,
                "Cn" : 9,
                "Ic" : 9,
                "Aa" : 9,
                "Cn" : 9,
                "Ic" : 9,
                "Bb" : 9,
                "An" : 9,
                "Ia" : 9,
                "Cc" : 9,
                "Bn" : 9,
                "Ib" : 9,
                "P2/m" : 10,
                "P2/m" : 10,
                "P2/m" : 10,
                "P21/m" : 11,
                "P21/m" : 11,
                "P21/m" : 11,
                "C2/m" : 12,
                "A2/m" : 12,
                "I2/m" : 12,
                "A2/m" : 12,
                "B2/m" : 12,
                "I2/m" : 12,
                "B2/m" : 12,
                "C2/m" : 12,
                "I2/m" : 12,
                "P2/c" : 13,
                "P2/n" : 13,
                "P2/a" : 13,
                "P2/a" : 13,
                "P2/n" : 13,
                "P2/b" : 13,
                "P2/b" : 13,
                "P2/n" : 13,
                "P2/c" : 13,
                "P21/c" : 14,
                "P21/n" : 14,
                "P21/a" : 14,
                "P21/a" : 14,
                "P21/n" : 14,
                "P21/b" : 14,
                "P21/b" : 14,
                "P21/n" : 14,
                "P21/c" : 14,
                "C2/c" : 15,
                "A2/n" : 15,
                "I2/a" : 15,
                "A2/a" : 15,
                "B2/n" : 15,
                "I2/b" : 15,
                "B2/b" : 15,
                "C2/n" : 15,
                "I2/c" : 15,
                "A2/a" : 15,
                "C2/n" : 15,
                "I2/c" : 15,
                "B2/b" : 15,
                "A2/n" : 15,
                "I2/a" : 15,
                "C2/c" : 15,
                "B2/n" : 15,
                "I2/b" : 15,
                "P222" : 16,
                "P2221" : 17,
                "P21212" : 18,
                "P212121" : 19,
                "C2221" : 20,
                "C222" : 21,
                "F222" : 22,
                "I222" : 23,
                "I212121" : 24,
                "Pmm2" : 25,
                "Pmc21" : 26,
                "Pcc2" : 27,
                "Pma2" : 28,
                "Pca21" : 29,
                "Pnc2" : 30,
                "Pmn21" : 31,
                "Pba2" : 32,
                "Pna21" : 33,
                "Pnn2" : 34,
                "Cmm2" : 35,
                "Cmc21" : 36,
                "Ccc2" : 37,
                "Amm2" : 38,
                "Abm2" : 39,
                "Ama2" : 40,
                "Aba2" : 41,
                "Fmm2" : 42,
                "Fdd2" : 43,
                "Imm2" : 44,
                "Iba2" : 45,
                "Ima2" : 46,
                "Pmmm" : 47,
                "Pnnn" : 48,
                "Pnnn" : 48,
                "Pccm" : 49,
                "Pban" : 50,
                "Pban" : 50,
                "Pmma" : 51,
                "Pnna" : 52,
                "Pmna" : 53,
                "Pcca" : 54,
                "Pbam" : 55,
                "Pccn" : 56,
                "Pbcm" : 57,
                "Pnnm" : 58,
                "Pmmn" : 59,
                "Pmmn" : 59,
                "Pbcn" : 60,
                "Pbca" : 61,
                "Pnma" : 62,
                "Cmcm" : 63,
                "Cmca" : 64,
                "Cmmm" : 65,
                "Cccm" : 66,
                "Cmma" : 67,
                "Ccca" : 68,
                "Ccca" : 68,
                "Fmmm" : 69,
                "Fddd" : 70,
                "Fddd" : 70,
                "Immm" : 71,
                "Ibam" : 72,
                "Ibca" : 73,
                "Imma" : 74,
                "P4" : 75,
                "P41" : 76,
                "P42" : 77,
                "P43" : 78,
                "I4" : 79,
                "I41" : 80,
                "P-4" : 81,
                "I-4" : 82,
                "P4/m" : 83,
                "P42/m" : 84,
                "P4/n" : 85,
                "P4/n" : 85,
                "P42/n" : 86,
                "P42/n" : 86,
                "I4/m" : 87,
                "I41/a" : 88,
                "I41/a" : 88,
                "P422" : 89,
                "P4212" : 90,
                "P4122" : 91,
                "P41212" : 92,
                "P4222" : 93,
                "P42212" : 94,
                "P4322" : 95,
                "P43212" : 96,
                "I422" : 97,
                "I4122" : 98,
                "P4mm" : 99,
                "P4bm" : 100,
                "P42cm" : 101,
                "P42nm" : 102,
                "P4cc" : 103,
                "P4nc" : 104,
                "P42mc" : 105,
                "P42bc" : 106,
                "I4mm" : 107,
                "I4cm" : 108,
                "I41md" : 109,
                "I41cd" : 110,
                "P-42m" : 111,
                "P-42c" : 112,
                "P-421m" : 113,
                "P-421c" : 114,
                "P-4m2" : 115,
                "P-4c2" : 116,
                "P-4b2" : 117,
                "P-4n2" : 118,
                "I-4m2" : 119,
                "I-4c2" : 120,
                "I-42m" : 121,
                "I-42d" : 122,
                "P4/mmm" : 123,
                "P4/mcc" : 124,
                "P4/nbm" : 125,
                "P4/nbm" : 125,
                "P4/nnc" : 126,
                "P4/nnc" : 126,
                "P4/mbm" : 127,
                "P4/mnc" : 128,
                "P4/nmm" : 129,
                "P4/nmm" : 129,
                "P4/ncc" : 130,
                "P4/ncc" : 130,
                "P42/mmc" : 131,
                "P42/mcm" : 132,
                "P42/nbc" : 133,
                "P42/nbc" : 133,
                "P42/nnm" : 134,
                "P42/nnm" : 134,
                "P42/mbc" : 135,
                "P42/mnm" : 136,
                "P42/nmc" : 137,
                "P42/nmc" : 137,
                "P42/ncm" : 138,
                "P42/ncm" : 138,
                "I4/mmm" : 139,
                "I4/mcm" : 140,
                "I41/amd" : 141,
                "I41/amd" : 141,
                "I41/acd" : 142,
                "I41/acd" : 142,
                "P3" : 143,
                "P31" : 144,
                "P32" : 145,
                "P-3" : 147,
                "P312" : 149,
                "P321" : 150,
                "P3112" : 151,
                "P3121" : 152,
                "P3212" : 153,
                "P3221" : 154,
                "P3m1" : 156,
                "P31m" : 157,
                "P3c1" : 158,
                "P31c" : 159,
                "P-31m" : 162,
                "P-31c" : 163,
                "P-3m1" : 164,
                "P-3c1" : 165,
                "P6" : 168,
                "P61" : 169,
                "P65" : 170,
                "P62" : 171,
                "P64" : 172,
                "P63" : 173,
                "P-6" : 174,
                "P6/m" : 175,
                "P63/m" : 176,
                "P622" : 177,
                "P6122" : 178,
                "P6522" : 179,
                "P6222" : 180,
                "P6422" : 181,
                "P6322" : 182,
                "P6mm" : 183,
                "P6cc" : 184,
                "P63cm" : 185,
                "P63mc" : 186,
                "P-6m2" : 187,
                "P-6c2" : 188,
                "P-62m" : 189,
                "P-62c" : 190,
                "P6/mmm" : 191,
                "P6/mcc" : 192,
                "P63/mcm" : 193,
                "P63/mmc" : 194,
                "P23" : 195,
                "F23" : 196,
                "I23" : 197,
                "P213" : 198,
                "I213" : 199,
                "Pm-3" : 200,
                "Pn-3" : 201,
                "Pn-3" : 201,
                "Fm-3" : 202,
                "Fd-3" : 203,
                "Fd-3" : 203,
                "Im-3" : 204,
                "Pa-3" : 205,
                "Ia-3" : 206,
                "P432" : 207,
                "P4232" : 208,
                "F432" : 209,
                "F4132" : 210,
                "I432" : 211,
                "P4332" : 212,
                "P4132" : 213,
                "I4132" : 214,
                "P-43m" : 215,
                "F-43m" : 216,
                "I-43m" : 217,
                "P-43n" : 218,
                "F-43c" : 219,
                "I-43d" : 220,
                "Pm-3m" : 221,
                "Pn-3n" : 222,
                "Pn-3n" : 222,
                "Pm-3n" : 223,
                "Pn-3m" : 224,
                "Pn-3m" : 224,
                "Fm-3m" : 225,
                "Fm-3c" : 226,
                "Fd-3m" : 227,
                "Fd-3m" : 227,
                "Fd-3c" : 228,
                "Fd-3c" : 228,
                "Im-3m" : 229,
                "Ia-3d" : 230,

                "P1":1,
                "P-1":2,
                "P121":3,
                "P1211":4,
                "C121":5,
                "P1m1":6,
                "P1c1":7,
                "C1m1":8,
                "C1c1":9,
                "P12/m1":10,
                "P121/m1":11,
                "C12/m1":12,
                "P12/c1":13,
                "P121/c1":14,
                "P121/a1":14,
                "C12/c1":15,
                "P222":16,
                "P2221":17,
                "P21212":18,
                "P212121":19,
                "C2221":20,
                "C222":21,
                "F222":22,
                "I222":23,
                "I212121":24,
                "Pmm2":25,
                "Pmc21":26,
                "Pcc2":27,
                "Pma2":28,
                "Pca21":29,
                "Pnc2":30,
                "Pmn21":31,
                "Pba2":32,
                "Pna21":33,
                "Pnn2":34,
                "Cmm2":35,
                "Cmc21":36,
                "Ccc2":37,
                "Amm2":38,
                "Abm2":39,
                "Aba2":41,
                "Fmm2":42,
                "Fdd2":43,
                "Imm2":44,
                "Iba2":45,
                "Ima2":46,
                "P2/m2/m2/m":47,
                "P2/n2/n2/n":48,
                "P2/b2/a2/n":50,
                "P21/m2/m2/a":51,
                "P2/n21/n2/a":52,
                "P2/m2/n21/a":53,
                "P21/c2/c2/a":54,
                "P21/b21/a2/m":55,
                "P21/c21/c2/n":56,
                "P2/b21/c21/m":57,
                "P21/n21/n2/m":58,
                "P21/m21/m2/n":59,
                "P21/b2/c21/n":60,
                "P21/b21/c21/a":61,
                "P21/n21/m21/a":62,
                "C2/m2/c21/m":63,
                "C2/m2/c21/a":64,
                "C2/m2/m2/m":65,
                "C2/c2/c2/m":66,
                "C2/m2/m2/a":67,
                "C2/c2/c2/a":68,
                "F2/m2/m2/m":69,
                "F2/d2/d2/d":70,
                "I2/m2/m2/m":71,
                "I2/b2/a2/m":72,
                "I21/b21/c21/a":73,
                "I21/m21/m21/a":74,
                "P4":75,
                "P41":76,
                "P42":77,
                "P43":78,
                "I4":79,
                "I41":80,
                "P-4":81,
                "I-4":82,
                "P4/m":83,
                "P42/m":84,
                "P4/n":85,
                "P42/n":86,
                "I4/m":87,
                "I41/a":88,
                "P4212":90,
                "P4122":91,
                "P41212":92,
                "P4322":95,
                "P43212":96,
                "I422":97,
                "I4122":98,
                "P4mm":99,
                "P4bm":100,
                "P42nm":102,
                "P4nc":104,
                "P42mc":105,
                "I4mm":107,
                "I4cm":108,
                "I41md":109,
                "I41cd":110,
                "P-42m":111,
                "P-42c":112,
                "P-421m":113,
                "P-421c":114,
                "P-4m2":115,
                "P-4c2":116,
                "P-4b2":117,
                "P-4n2":118,
                "I-4m2":119,
                "I-4c2":120,
                "I-42m":121,
                "I-42d":122,
                "P4/m2/m2/m":123,
                "P4/m2/c2/c":124,
                "P4/n2/b2/m":125,
                "P4/n2/n2/c":126,
                "P4/m21/b2/m":127,
                "P4/m21/n2/c":128,
                "P4/n21/m2/m":129,
                "P4/n21/c2/c":130,
                "P42/m2/m2/c":131,
                "P42/m2/c2/m":132,
                "P42/n2/b2/c":133,
                "P42/n2/n2/m":134,
                "P42/m21/b2/c":135,
                "P42/m21/n2/m":136,
                "P42/n21/m2/c":137,
                "P42/n21/c2/m":138,
                "I4/m2/m2/m":139,
                "I4/m2/c2/m":140,
                "I41/a2/m2/d":141,
                "I41/a2/c2/d":142,
                "P3":143,
                "P31":144,
                "P32":145,
                "R3":146,
                "P-3":147,
                "R-3":148,
                "P312":149,
                "P321":150,
                "P3112":151,
                "P3121":152,
                "P3212":153,
                "P3221":154,
                "R32":155,
                "P3m1":156,
                "P31m":157,
                "P3c1":158,
                "P31c":159,
                "R3m":160,
                "R3c":161,
                "P-312/m":162,
                "P-312/c":163,
                "P-32/m1":164,
                "P-32/c1":165,
                "R-32/m":166,
                "R-32/c":167,
                "P61":169,
                "P65":170,
                "P63":173,
                "P-6":174,
                "P6/m":175,
                "P63/m":176,
                "P6122":178,
                "P6522":179,
                "P6222":180,
                "P6422":181,
                "P6322":182,
                "P6mm":183,
                "P63cm":185,
                "P63mc":186,
                "P-6m2":187,
                "P-6c2":188,
                "P-62m":189,
                "P-62c":190,
                "P6/m2/m2/m":191,
                "P6/m2/c2/c":192,
                "P63/m2/c2/m":193,
                "P63/m2/m2/c":194,
                "P23":195,
                "F23":196,
                "I23":197,
                "P213":198,
                "I213":199,
                "P2/m-3":200,
                "P2/n-3":201,
                "F2/m-3":202,
                "F2/d-3":203,
                "I2/m-3":204,
                "P21/a-3":205,
                "I21/a-3":206,
                "P4232":208,
                "F4132":210,
                "I432":211,
                "P4332":212,
                "P4132":213,
                "I4132":214,
                "P-43m":215,
                "F-43m":216,
                "I-43m":217,
                "P-43n":218,
                "F-43c":219,
                "I-43d":220,
                "P4/m-32/m":221,
                "P4/n-32/n":222,
                "P42/m-32/n":223,
                "P42/n-32/m":224,
                "F4/m-32/m":225,
                "F4/m-32/c":226,
                "F41/d-32/m":227,
                "I4/m-32/m":229,
                "I41/a-32/d":230,}
        try:
            return HMdict[HM_str.replace(" ","")]
        except Exception as e: 
#            print e
            raise SystemExit
#            return 1




def reduce_atoms(all_eq_pos,labels,cell=np.identity(3,dtype=float),thresh=0.1,index=False):
    '''shift all atoms back into cell if need be'''
    all_eq_pos%=1.0
    '''make a copy to keep in crystal coords'''
    all_eq_pos_crys = np.copy(all_eq_pos)
    '''distances in alat'''
    all_eq_pos = (cell.dot(all_eq_pos.T)).T
    '''distance between each atom pair'''
    dist = periodic_dist_func(all_eq_pos,all_eq_pos,cell=cell)
    '''mask for duplicate positions'''
    mask = np.ones(all_eq_pos.shape[0],dtype=bool)
    '''for self mask'''
    xr = list(range(all_eq_pos.shape[0]))
    idx = np.array(xr,dtype=int)
    '''loop over each atomic position'''
    for i in xr:
        '''if already mask don't mask on this indice'''
        if not mask[i]:
            continue
        '''self mask'''
        self_bool = idx==i
        '''distance mask for duplicates'''
        dist_bool = dist[i]>thresh
        '''if distance is less than the threshold and it's the same atom or if it's'''
        '''a different atom and the position is greater than the threshold allow'''
        '''all masks combined'''
        mask  *= np.logical_or(dist_bool,self_bool)

    '''reduce positions to non duplicates..back to fractional coords'''
    all_eq_pos = all_eq_pos_crys[mask]
    '''reduce labels'''
    labels_arr=np.array(labels)[mask]

    '''sort by species label'''
    spec_sort = np.argsort(labels_arr)
    labels_arr=labels_arr[spec_sort]
    all_eq_pos=all_eq_pos[spec_sort]

    if index==True:
        return mask
    
    else: return all_eq_pos,labels_arr

def periodic_dist_func(X,Y,cell=np.identity(3,dtype=float)):
    '''
    needed to check for distance atoms of a system 
    with periodic boundary conditions
    '''


    try:
        Xt=np.copy(Y)
        Xtc=np.copy(Y)
        dist=np.zeros((X.shape[0],Xt.shape[0],8))

        dist[:,:,0] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,0.0,0.0],]).T)).T

        dist[:,:,1] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[0.0,1.0,0.0],]).T)).T

        dist[:,:,2] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[0.0,0.0,1.0],]).T)).T

        dist[:,:,3] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,0.0,1.0],]).T)).T

        dist[:,:,4] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,1.0,0.0],]).T)).T

        dist[:,:,5] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[0.0,1.0,1.0],]).T)).T

        dist[:,:,6] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,1.0,1.0],]).T)).T

        dist[:,:,7] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')

        return np.amin(dist,axis=2)

    except Exception as e:
        print(e)
        




