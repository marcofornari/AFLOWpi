import AFLOWpi
import numpy as np
import os
import subprocess
import StringIO
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
        self.sg_num    = self.__get_sg_num()
        self.ibrav     = self.__ibrav_from_sg_number()
        self.iso_pr_car= ''
        
        self.iso_basis = self.__get_iso_basis()
        self.iso_pr_car= self.__get_iso_cart()

        input_dict = AFLOWpi.retr._splitInput(self.input)
        print input_dict["&system"]["nat"]

    def qe_output(self,input_str,accuracy=0.001):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()


        self.input     =  AFLOWpi.prep._removeComments(input_str)
        self.accuracy  = accuracy

        self.output    = self.__get_isotropy_output(qe_output=True)
        self.cif       = self.qe2cif()
        self.sg_num    = self.__get_sg_num()
        self.ibrav     = self.__ibrav_from_sg_number()
        self.iso_pr_car= ''

        self.iso_basis = self.__get_iso_basis()
        self.iso_pr_car= self.__get_iso_cart()

    def cif_input(self,input_str):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()

        self.cif       = input_str
        self.output    = input_str
        self.input     = self.__skeleton_input()
        self.accuracy  = 0.0
        self.sg_num    = self.__get_sg_num()
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
        return re.findall('# CIF file.*\n(?:.*\n)*',self.output)[0]

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
                    positions.append(' '.join(split_pos[1:4]))
                except:
                    pass
            positions='\n'.join(positions)

        else:
            cell_matrix = AFLOWpi.retr.getCellMatrixFromInput(self.input,string=False)*0.529177249
            print cell_matrix

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
        centering='P'

        in_str = self.__generate_isotropy_input_from_qe_data(output=qe_output)
        ISODATA = os.path.join(AFLOWpi.__path__[0],'ISOTROPY')
        os.putenv('ISODATA',ISODATA+'/') 
        findsym_path = os.path.join(ISODATA,'findsym')


        try:
            find_sym_process = subprocess.Popen(findsym_path,stdin=subprocess.PIPE,stdout=subprocess.PIPE,)
            output = find_sym_process.communicate(input=in_str)[0]
            self.output=output
        except:
            print find_sym_process.returncode

        return output

    def __get_sg_num(self):
        try:
            sg_info = re.findall('Space Group\s*([0-9]*)\s*([\w-]*)\s*([\w-]*)',self.output)[0]
            self.sgn = int(sg_info[0])
        except:
            try:
                sg_info = re.findall('_symmetry_Int_Tables_number\s*(\d+)',self.cif)[0]
                self.sgn = int(sg_info)
            except:
                try:
                    sg_nam = re.compile('''_symmetry_space_group_name_H-M\s*['"](.*)['"]''')
                    self.sgn = self.__HM2Num(sg_nam.findall(self.cif)[0])
                    return self.sgn

                except Exception,e:
                    print e
#                    print self.cif
#                    raise SystemExit

    def __get_iso_cart(self):
            search = 'Lattice vectors in cartesian coordinates:\s*\n(.*\n.*\n.*)\n'
            std_prim_basis_str = re.findall(search,self.output)[0]
            self.iso_pr_car    = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)#*

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

        modifier = AFLOWpi.qe.regex.cell_parameters(self.input,return_which='modifier',).lower()

        origin = self.__conv_from_cif()

        origin = np.asarray([float(i) for i in origin])
        self.origin=origin

        std_prim_basis_str = re.findall('Vectors a,b,c:\s*\n(.*\n.*\n.*)\n',self.output)[0]
        self.iso_conv = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)

        input_dict = AFLOWpi.retr._splitInput(self.input)
        if 'CELL_PARAMETERS' not in input_dict:
            prim_in = AFLOWpi.retr.getCellMatrixFromInput(self.input)
            try:
                prim_in=AFLOWpi.retr._cellStringToMatrix(prim_in)
            except:
                pass
        else:
            if input_dict['CELL_PARAMETERS']["__content__"]=='':
                try:
                    prim_in = AFLOWpi.retr.getCellMatrixFromInput(self.input)
                except:
                    pass
            try:
                prim_in=AFLOWpi.retr._cellStringToMatrix(input_dict['CELL_PARAMETERS']['__content__'])
            except Exception,e:
                print e
                raise SystemExit

        self.iso_basis=prim_in

        a=np.array([[self.conv_a,],
                       [self.conv_b,],
                       [self.conv_c,],])


        a=np.abs((self.iso_conv).dot(self.iso_basis))

        if self.ibrav in [8,9,10,11]:
            self.conv_a=np.sum(a[:,0])
            self.conv_b=np.sum(a[:,1])
            self.conv_c=np.sum(a[:,2])


        prim_abc = re.findall('Lattice parameters, a,b,c,alpha,beta,gamma.*\n(.*)\n',self.output)[0].split()
        prim_abc=map(float,prim_abc)

        return self.iso_basis


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
            print 'SG num not found:',self.sgn
            raise SystemExit



    def convert(self,ibrav=True,thresh=0.01):
        inputString=''

        for thresh in [0.001,0.01,0.1]:
            inputString_new= self.cif2qe(thresh=thresh)
            isotropy_chk    = AFLOWpi.prep.isotropy()
            isotropy_chk.qe_input(inputString_new)


            inputString=inputString_new
            break
            if self.sgn==isotropy_chk.sgn:

                break
            else:
                print self.sgn,isotropy_chk.sgn
        return inputString


    def cif2qe(self,thresh=0.01):

        input_dict = AFLOWpi.retr._splitInput(self.input)
        '''grab the conventional -> primitive conversion matrix for this ibrav'''
        if self.ibrav in [4,5]:
            t_ibrav=1
        else:
            t_ibrav=self.ibrav



        '''conventional to qe convention primitive lattice vec'''
        convert = AFLOWpi.retr.abc2free(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,
                                        gamma=90.0,ibrav=t_ibrav,returnString=False)

        if self.ibrav==5:

             in_hex = AFLOWpi.retr.abc2free(a=self.conv_a,b=1.0,c=self.conv_c,alpha=90.0,
                                            beta=90.0,gamma=120.0,ibrav=4,
                                            returnString=False)
            # print in_hex
             beta=np.sqrt(3.0+(self.conv_c/self.conv_a)**2.0)
             self.conv_a = self.conv_a*beta/3.0
             self.conv_b = self.conv_a
             self.conv_c = self.conv_a
             self.conv_alpha =  2.0*np.arcsin((3.0/(2.0*beta)))*(180.0/np.pi)
             self.conv_beta  = self.conv_alpha
             self.conv_gamma = self.conv_alpha
             in_rho = AFLOWpi.retr.abc2free(a=self.conv_a,b=self.conv_b,c=self.conv_c,
                                            alpha=self.conv_alpha,beta=self.conv_alpha,
                                            gamma=self.conv_gamma,ibrav=5,
                                            returnString=False)

             convert = in_hex.dot(np.linalg.inv(in_rho))
             '''hex to rho'''
             convert = np.array([[-1.,  1., -0.,],
                                 [ 1.,  0., -1.,],
                                 [ 1.,  1.,  1.,],])
             
             convert=np.linalg.inv(np.around(convert,decimals=1))

                                    

        ins= self.cif.lower()
        '''grab the symmetry operations from the text in the cif'''
        re_symops=re.compile(r'_space_group_symop_operation_xyz\s*\n((?:\s*\d+.*\n)+)\s*',re.M)

        symops = [x for x in re_symops.findall(ins)[0].split('\n') if (len(x.strip())!=0 )]

        re_sym_ops_remove_numbering = re.compile('\d+\s+(.*)')
        for i in  range(len(symops)):
            symops[i] = re_sym_ops_remove_numbering.findall(symops[i])[0]

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

        spec_lab =  loop_list.index('_atom_site_label')

        x_loc =  loop_list.index('_atom_site_fract_x')
        y_loc =  loop_list.index('_atom_site_fract_y')
        z_loc =  loop_list.index('_atom_site_fract_z')

        '''get positions from input'''
        positions = [map(str.strip,x.split()) for x in  atom_pos[1].split('\n') if len(x.strip())!=0]

        pos_array=np.zeros((len(positions),3))
        
        labels=[]
        for i in range(len(positions)):
            pos_array[i][0] = float(positions[i][x_loc])
            pos_array[i][1] = float(positions[i][y_loc])
            pos_array[i][2] = float(positions[i][z_loc])

            labels.append(positions[i][spec_lab].strip('0123456789').title())

        '''shift positions to between 0.0 and 1.0'''
        pos_array%=1.0

        '''if positions is 1/3 or 2/3 use more precision on position'''
        pos_array[np.where(np.isclose(pos_array-(1.0/3.0),0.0,rtol=1e-04, atol=1e-05))]=1.0/3.0
        pos_array[np.where(np.isclose(pos_array-(2.0/3.0),0.0,rtol=1e-04, atol=1e-05))]=2.0/3.0

        '''there will by natm*numSymOps'''
        labels = labels*len(operations)
        all_eq_pos=np.zeros((pos_array.shape[0]*operations.shape[0],3))

        '''duplicate atoms by symmetry operations'''
        for i in xrange(operations.shape[0]):

            temp_pos=np.copy(pos_array)
            temp_pos   = temp_pos.dot(operations[i])
            '''translation symmetry'''
            temp_pos[:,0]+=shift[i][0]
            temp_pos[:,1]+=shift[i][1]
            temp_pos[:,2]+=shift[i][2]

            '''add atoms generated by symmetry operations to the list of atoms'''
            all_eq_pos[i*pos_array.shape[0]:(i+1)*pos_array.shape[0]]=temp_pos


#        '''shift origin to zero'''
#        all_eq_pos-=self.origin

        if self.sgn in [38,39,40,41]:
            all_eq_pos=all_eq_pos[:,[2,1,0]]
            all_eq_pos[:,0]*=-1.0
            temp_c = self.conv_c
            self.conv_c=self.conv_a
            self.conv_a=temp_c



            convert = AFLOWpi.retr.abc2free(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,
                                            gamma=self.conv_gamma,ibrav=self.ibrav,returnString=False)


        '''make unique a or b monoclinic into unique c'''
        if self.ibrav in [12,13]:
            if np.isclose(self.conv_beta,self.conv_gamma):
                #all_eq_pos=all_eq_pos[:,[1,2,0]]
                temp_a=self.conv_b
                temp_b=self.conv_c
                temp_c=self.conv_a
                temp_alpha=self.conv_gamma
                temp_beta =self.conv_alpha
                temp_gamma=self.conv_beta
                self.conv_a=temp_a
                self.conv_b=temp_b
                self.conv_c=temp_c
                self.conv_alpha = temp_alpha
                self.conv_beta  = temp_beta
                self.conv_gamma = temp_gamma

            elif np.isclose(self.conv_alpha,self.conv_gamma):
                all_eq_pos=all_eq_pos[:,[2,0,1]]

                print self.conv_alpha
                print self.conv_beta                    
                print self.conv_gamma
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
        '''######################'''
        '''remove duplicate atoms'''
        '''######################'''

        if self.ibrav==13:
            self.ibrav=12
            two_one = np.cos(self.conv_gamma/180.0*np.pi)
            two_two = np.sin(self.conv_gamma/180.0*np.pi)

            convert=np.matrix((
                    ( 0.5, 0.0,0.5),
                    ( 0.5, 0.0,-0.5),
                    ( two_one,two_two,0.0),
#                    ( 0.5, 0.0, 0.5,),
                               )).T
            
        if self.ibrav==12:
            convert=np.matrix(((1.0, 0.0,0.0,),
#                               (1.0*np.cos(self.conv_gamma), 1.0*np.sin(self.conv_gamma), 0.0,),
                               (0.0, 1.0, 0.0,),
                               (0.0, 0.0, 1.0,),))



        '''shift all atoms back into cell if need be'''
#        all_eq_pos%=1.0
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
            print self.conv_a
            print self.conv_b
            print self.conv_c
            print self.conv_alpha
            print self.conv_beta                    
            print self.conv_gamma

            two_one = np.cos(self.conv_gamma/180.0*np.pi)*self.conv_b
            two_two = np.sin(self.conv_gamma/180.0*np.pi)*self.conv_b
#            to_conv=np.matrix((
#                    (0.5*self.conv_a, 0.0        ,0.5*self.conv_c,),
#                    (1.0*two_one,1.0*two_two, 0.0            ,),
#                    (-0.5*self.conv_a, 0.0        ,-0.5*self.conv_c,),))
            to_conv=np.matrix(((1.0*self.conv_a, 0.0,0.0,),
                               (0.0, 1.0*self.conv_b, 0.0,),
                               (0.0, 0.0, 1.0*self.conv_c,),))


            # to_conv=np.matrix(((1.0*self.conv_a ,0.0,0.0,),
            #                    (1.0*np.cos(self.conv_gamma)*self.conv_b, 1.0*np.sin(self.conv_gamma)*self.conv_b, 0.0,),
            #                      (0.0,0.0,1.0*self.conv_c,),)).T
        if self.ibrav==12:
            to_conv=np.matrix(((1.0*self.conv_a, 0.0,0.0,),
                               (0.0, 1.0*self.conv_b, 0.0,),
                               (0.0, 0.0, 1.0*self.conv_c,),))



        '''reduce equivilent atoms'''
        all_eq_pos,labels_arr = reduce_atoms(all_eq_pos,labels,cell=to_conv,thresh=thresh)


        '''form atomic positions card for QE input'''
        atm_pos_str=""
        for i in xrange(all_eq_pos.shape[0]):
            atm_pos_str+= ('%3.3s'% labels_arr[i]) +(' % 9.9f % 9.9f % 9.9f '%tuple(all_eq_pos[i].tolist()))+"\n"





        '''assign values of cell to A,B,C,cosAB,cosAC,cosBC in QE input file'''
        input_dict['&system']['ibrav']=self.ibrav

        input_dict['&system']['A']=self.conv_a
        input_dict['&system']['B']=self.conv_b
        input_dict['&system']['C']=self.conv_c

        input_dict['&system']['cosAB']=np.cos(self.conv_gamma/180.0*np.pi)
        input_dict['&system']['cosAC']=np.cos(self.conv_beta/180.0*np.pi)
        input_dict['&system']['cosBC']=np.cos(self.conv_alpha/180.0*np.pi)


        
                
        try:
            del input_dict['CELL_PARAMETERS']
        except:
            pass

    
        
        input_dict['&system']['nat']=all_eq_pos.shape[0]
        
        '''add newly transformed atomic positions to QE convention input'''
        try:
            input_dict['ATOMIC_POSITIONS']['__content__']=atm_pos_str        
        except:
            input_dict['ATOMIC_POSITIONS']={}
            input_dict['ATOMIC_POSITIONS']["__modifier__"]="{crystal}"
            input_dict['ATOMIC_POSITIONS']['__content__']=atm_pos_str        

        qe_convention_input=AFLOWpi.retr._joinInput(input_dict)
        '''convert to celldm'''
        qe_convention_input = AFLOWpi.prep._transformInput(qe_convention_input)


        return qe_convention_input


    def __HM2Num(self,HM_str):
       HM_dict={"P 1":1,
                "P -1":2,
                "P 1 2 1":3,
                "P 1 21 1":4,
                "C 1 2 1":5,
                "P 1 m 1":6,
                "P 1 c 1":7,
                "C 1 m 1":8,
                "C 1 c 1":9,
                "P 1 2/m 1":10,
                "P 1 21/m 1":11,
                "C 1 2/m 1":12,
                "P 1 2/c 1":13,
                "P 1 21/c 1":14,
                "P 1 21/a 1":14,
                "C 1 2/c 1":15,
                "P 2 2 2":16,
                "P 2 2 21":17,
                "P 21 21 2":18,
                "P 21 21 21":19,
                "C 2 2 21":20,
                "C 2 2 2":21,
                "F 2 2 2":22,
                "I 2 2 2":23,
                "I 21 21 21":24,
                "P m m 2":25,
                "P m c 21":26,
                "P c c 2":27,
                "P m a 2":28,
                "P c a 21":29,
                "P n c 2":30,
                "P m n 21":31,
                "P b a 2":32,
                "P n a 21":33,
                "P n n 2":34,
                "C m m 2":35,
                "C m c 21":36,
                "C c c 2":37,
                "A m m 2":38,
                "A b m 2":39,
                "A b a 2":41,
                "F m m 2":42,
                "F d d 2":43,
                "I m m 2":44,
                "I b a 2":45,
                "I m a 2":46,
                "P 2/m 2/m 2/m":47,
                "P 2/n 2/n 2/n (origin choice 2)":48,
                "P 2/b 2/a 2/n (origin choice 2)":50,
                "P 21/m 2/m 2/a":51,
                "P 2/n 21/n 2/a":52,
                "P 2/m 2/n 21/a":53,
                "P 21/c 2/c 2/a":54,
                "P 21/b 21/a 2/m":55,
                "P 21/c 21/c 2/n":56,
                "P 2/b 21/c 21/m":57,
                "P 21/n 21/n 2/m":58,
                "P 21/m 21/m 2/n (origin choice 2)":59,
                "P 21/b 2/c 21/n":60,
                "P 21/b 21/c 21/a":61,
                "P 21/n 21/m 21/a":62,
                "C 2/m 2/c 21/m":63,
                "C 2/m 2/c 21/a":64,
                "C 2/m 2/m 2/m":65,
                "C 2/c 2/c 2/m":66,
                "C 2/m 2/m 2/a":67,
                "C 2/c 2/c 2/a (origin choice 2)":68,
                "F 2/m 2/m 2/m":69,
                "F 2/d 2/d 2/d (origin choice 2)":70,
                "I 2/m 2/m 2/m":71,
                "I 2/b 2/a 2/m":72,
                "I 21/b 21/c 21/a":73,
                "I 21/m 21/m 21/a":74,
                "P 4":75,
                "P 41":76,
                "P 42":77,
                "P 43":78,
                "I 4":79,
                "I 41":80,
                "P -4":81,
                "I -4":82,
                "P 4/m":83,
                "P 42/m":84,
                "P 4/n (origin choice 2)":85,
                "P 42/n (origin choice 2)":86,
                "I 4/m":87,
                "I 41/a (origin choice 2)":88,
                "P 4 21 2":90,
                "P 41 2 2":91,
                "P 41 21 2":92,
                "P 43 2 2":95,
                "P 43 21 2":96,
                "I 4 2 2":97,
                "I 41 2 2":98,
                "P 4 m m":99,
                "P 4 b m":100,
                "P 42 n m":102,
                "P 4 n c":104,
                "P 42 m c":105,
                "I 4 m m":107,
                "I 4 c m":108,
                "I 41 m d":109,
                "I 41 c d":110,
                "P -4 2 m":111,
                "P -4 2 c":112,
                "P -4 21 m":113,
                "P -4 21 c":114,
                "P -4 m 2":115,
                "P -4 c 2":116,
                "P -4 b 2":117,
                "P -4 n 2":118,
                "I -4 m 2":119,
                "I -4 c 2":120,
                "I -4 2 m":121,
                "I -4 2 d":122,
                "P 4/m 2/m 2/m":123,
                "P 4/m 2/c 2/c":124,
                "P 4/n 2/b 2/m (origin choice 2)":125,
                "P 4/n 2/n 2/c (origin choice 2)":126,
                "P 4/m 21/b 2/m":127,
                "P 4/m 21/n 2/c":128,
                "P 4/n 21/m 2/m (origin choice 2)":129,
                "P 4/n 21/c 2/c (origin choice 2)":130,
                "P 42/m 2/m 2/c":131,
                "P 42/m 2/c 2/m":132,
                "P 42/n 2/b 2/c (origin choice 2)":133,
                "P 42/n 2/n 2/m (origin choice 2)":134,
                "P 42/m 21/b 2/c":135,
                "P 42/m 21/n 2/m":136,
                "P 42/n 21/m 2/c (origin choice 2)":137,
                "P 42/n 21/c 2/m (origin choice 2)":138,
                "I 4/m 2/m 2/m":139,
                "I 4/m 2/c 2/m":140,
                "I 41/a 2/m 2/d (origin choice 2)":141,
                "I 41/a 2/c 2/d (origin choice 2)":142,
                "P 3":143,
                "P 31":144,
                "P 32":145,
                "R 3 (hexagonal axes)":146,
                "P -3":147,
                "R -3 (hexagonal axes)":148,
                "P 3 1 2":149,
                "P 3 2 1":150,
                "P 31 1 2":151,
                "P 31 2 1":152,
                "P 32 1 2":153,
                "P 32 2 1":154,
                "R 3 2 (hexagonal axes)":155,
                "P 3 m 1":156,
                "P 3 1 m":157,
                "P 3 c 1":158,
                "P 3 1 c":159,
                "R 3 m (hexagonal axes)":160,
                "R 3 c (hexagonal axes)":161,
                "P -3 1 2/m":162,
                "P -3 1 2/c":163,
                "P -3 2/m 1":164,
                "P -3 2/c 1":165,
                "R -3 2/m (hexagonal axes)":166,
                "R -3 2/c (hexagonal axes)":167,
                "P 61":169,
                "P 65":170,
                "P 63":173,
                "P -6":174,
                "P 6/m":175,
                "P 63/m":176,
                "P 61 2 2":178,
                "P 65 2 2":179,
                "P 62 2 2":180,
                "P 64 2 2":181,
                "P 63 2 2":182,
                "P 6 m m":183,
                "P 63 c m":185,
                "P 63 m c":186,
                "P -6 m 2":187,
                "P -6 c 2":188,
                "P -6 2 m":189,
                "P -6 2 c":190,
                "P 6/m 2/m 2/m":191,
                "P 6/m 2/c 2/c":192,
                "P 63/m 2/c 2/m":193,
                "P 63/m 2/m 2/c":194,
                "P 2 3":195,
                "F 2 3":196,
                "I 2 3":197,
                "P 21 3":198,
                "I 21 3":199,
                "P 2/m -3":200,
                "P 2/n -3 (origin choice 2)":201,
                "F 2/m -3":202,
                "F 2/d -3 (origin choice 2)":203,
                "I 2/m -3":204,
                "P 21/a -3":205,
                "I 21/a -3":206,
                "P 42 3 2":208,
                "F 41 3 2":210,
                "I 4 3 2":211,
                "P 43 3 2":212,
                "P 41 3 2":213,
                "I 41 3 2":214,
                "P -4 3 m":215,
                "F -4 3 m":216,
                "I -4 3 m":217,
                "P -4 3 n":218,
                "F -4 3 c":219,
                "I -4 3 d":220,
                "P 4/m -3 2/m":221,
                "P 4/n -3 2/n (origin choice 2)":222,
                "P 42/m -3 2/n":223,
                "P 42/n -3 2/m (origin choice 2)":224,
                "F 4/m -3 2/m":225,
                "F 4/m -3 2/c":226,
                "F 41/d -3 2/m (origin choice 2)":227,
                "I 4/m -3 2/m":229,
                "I 41/a -3 2/d":230,}

       return HM_dict[HM_str]




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
    xr = xrange(all_eq_pos.shape[0])
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

    print 
    try:
        Xt=np.copy(Y)
        Xtc=np.copy(Y)
        dist=np.zeros((X.shape[0],Xt.shape[0],8))

        dist[:,:,0] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,0.0,0.0],]).T)).T
#                  Xt=Xtc-np.array([1.0,0.0,0.0]).dot(cell)
        dist[:,:,1] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[0.0,1.0,0.0],]).T)).T
#                  Xt=Xtc-np.array([0.0,1.0,0.0]).dot(cell)
        dist[:,:,2] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[0.0,0.0,1.0],]).T)).T
#                  Xt=Xtc-np.array([0.0,0.0,1.0]).dot(cell)
        dist[:,:,3] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,0.0,1.0],]).T)).T
#                  Xt=Xtc-np.array([1.0,0.0,1.0]).dot(cell)
        dist[:,:,4] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,1.0,0.0],]).T)).T
#                  Xt=Xtc-np.array([1.0,1.0,0.0]).dot(cell)            
        dist[:,:,5] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[0.0,1.0,1.0],]).T)).T
#                  Xt=Xtc-np.array([0.0,1.0,1.0]).dot(cell)            
        dist[:,:,6] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')
        Xt=Xtc-(cell.dot(np.array([[1.0,1.0,1.0],]).T)).T
#                  Xt=Xtc-np.array([1.0,1.0,1.0]).dot(cell)
        dist[:,:,7] = scipy.spatial.distance.cdist(X, Xt, 'euclidean')

        return np.amin(dist,axis=2)

    except Exception,e:
        print e
        raise SystemExit


