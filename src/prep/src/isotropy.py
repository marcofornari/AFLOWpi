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



    def qe_input(self,input_str,accuracy=0.001):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()


        self.input     =  AFLOWpi.prep._removeComments(input_str)
        self.accuracy  = accuracy

        self.output    = self.__get_isotropy_output(qe_output=False)
        self.sg_num    = self.__get_sg_num()
        self.ibrav     = self.__ibrav_from_sg_number()
        self.cif       = self.qe2cif()
        self.iso_pr_car= ''

        
        self.iso_basis = self.__get_iso_basis()
        self.iso_pr_car= self.__get_iso_cart()


    def qe_output(self,input_str,accuracy=0.001):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()


        self.input     =  AFLOWpi.prep._removeComments(input_str)
        self.accuracy  = accuracy

        self.output    = self.__get_isotropy_output(qe_output=True)
        self.sg_num    = self.__get_sg_num()
        self.ibrav     = self.__ibrav_from_sg_number()
        self.cif       = self.qe2cif()
        self.iso_pr_car= ''

        
        self.iso_basis = self.__get_iso_basis()
        self.iso_pr_car= self.__get_iso_cart()

    def cif_input(self,input_str):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()


        self.input     = self.__skeleton_input()
        self.cif       = input_str
        self.accuracy  = 0.0

        self.output    = input_str
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
            positions = AFLOWpi.retr._getPositions(self.input,matrix=False)
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
            sg_info = re.findall('_symmetry_Int_Tables_number\s*(\d+)',self.cif)[0]
            self.sgn = int(sg_info)


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
            return 5

        elif self.sgn in [146,148,155,160,161,166,167]:
            return 4
                     
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



    def convert(self,ibrav=True):
        return self.cif2qe()


    def cif2qe(self):

        input_dict = AFLOWpi.retr._splitInput(self.input)
        '''grab the conventional -> primitive conversion matrix for this ibrav'''
        convert = AFLOWpi.retr.abc2free(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,gamma=90.0,ibrav=self.ibrav,returnString=False)

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
         
        '''make unique a or b monoclinic into unique c'''
        if self.ibrav in [12,13]:
            if np.isclose(self.conv_beta,self.conv_gamma):
                all_eq_pos=all_eq_pos[:,[1,2,0]]
            if np.isclose(self.conv_alpha,self.conv_gamma):
                all_eq_pos=all_eq_pos[:,[2,0,1]]


        '''######################'''
        '''remove duplicate atoms'''
        '''######################'''
        '''convert xyz to crystal coords'''
        try:
            all_eq_pos=(np.linalg.inv(convert.getA()).T.dot(all_eq_pos.T)).T
        except:
            all_eq_pos=(np.linalg.inv(convert).T.dot(all_eq_pos.T)).T

        '''shift all atoms back into cell if need be'''
        all_eq_pos%=1.0
        '''distance between each atom pair'''
        dist = scipy.spatial.distance.cdist(all_eq_pos,all_eq_pos)
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
            dist_bool = dist[i]>0.001
            '''if distance is less than the threshold and it's the same atom or if it's'''
            '''a different atom and the position is greater than the threshold allow'''
            '''all masks combined'''
            mask  *= np.logical_or(dist_bool,self_bool)

        '''reduce positions to non duplicates'''
        all_eq_pos = all_eq_pos[mask]
        '''reduce labels'''
        labels_arr=np.array(labels)[mask]

        '''sort by species label'''
        spec_sort = np.argsort(labels_arr)
        labels_arr=labels_arr[spec_sort]
        all_eq_pos=all_eq_pos[spec_sort]

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

        '''make unique a or b monoclinic into unique c'''
        if self.ibrav in [12,13]:
            if np.isclose(self.conv_alpha,self.conv_gamma):
                input_dict['&system']['A']=self.conv_b
                input_dict['&system']['B']=self.conv_c
                input_dict['&system']['C']=self.conv_a
                input_dict['&system']['cosBC']=np.cos(self.conv_gamma/180.0*np.pi)
                input_dict['&system']['cosAC']=np.cos(self.conv_alpha/180.0*np.pi)
                input_dict['&system']['cosAB']=np.cos(self.conv_beta/180.0*np.pi)
            if np.isclose(self.conv_beta,self.conv_gamma):
                input_dict['&system']['A']=self.conv_c
                input_dict['&system']['B']=self.conv_a
                input_dict['&system']['C']=self.conv_b
                input_dict['&system']['cosBC']=np.cos(self.conv_beta/180.0*np.pi)
                input_dict['&system']['cosAC']=np.cos(self.conv_gamma/180.0*np.pi)
                input_dict['&system']['cosAB']=np.cos(self.conv_alpha/180.0*np.pi)
        
                
        try:
            del input_dict['CELL_PARAMETERS']
        except:
            pass

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
