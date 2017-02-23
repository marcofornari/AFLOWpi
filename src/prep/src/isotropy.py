import AFLOWpi
import numpy
import os
import subprocess
import StringIO
import re
import copy



numpy.set_printoptions(precision=4, threshold=200, edgeitems=200, linewidth=250, suppress=True)                   
class isotropy():

    def __init__(self,input_str,accuracy=0.001,output=False):

        if os.path.exists(input_str):
            with open(input_str,'r') as fo:
                input_str = fo.read()
                

        self.input     =  AFLOWpi.prep._removeComments(input_str)
        self.accuracy  = accuracy
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
        self.axes_flip = ''
        self.output    = self.__get_isotropy_output(qe_output=output)
        self.sg_num    = self.__get_sg_num()
        self.ibrav     = self.__ibrav_from_sg_number()
        self.cif       = self.get_cif()
#        self.iso_basis = ''
        self.iso_pr_car= ''

        self.conv      = ''
#        if output==False:
        try:
            self.iso_basis = self.__get_iso_basis()

            self.iso_pr_car= self.__get_iso_cart()
            self.conv      = self.__convert_to_ibrav_matrix()
            

#            self.qe   = self.convert()
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            pass
#        self.a
#        self.b
#        self.c
#        self.iso_co_car=

#        raise SystemExit


    def get_cif(self):
        return re.findall('# CIF file.*\n(?:.*\n)*',self.output)[0]

    def __generate_isotropy_input_from_qe_data(self,output=False):
        if output:
#            cell_matrix = AFLOWpi.retr.getCellMatrixFromOutput(self.input,string=False)*0.529177249
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
#            print positions
#            print cm_string
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
#        print labels
        in_list={}
        for centering in ['P',]:
            isotropy_input_str='input file for isotropy generated from pwscf input by AFLOWpi\n'    
            isotropy_input_str+='%s\n'%self.accuracy
            isotropy_input_str+='1\n'
            isotropy_input_str+=cm_string
            isotropy_input_str+='2\n'
            isotropy_input_str+=centering+'\n'
            isotropy_input_str+='%s\n'%num_atoms
            isotropy_input_str+=' '.join([str(j) for j in labels])+'\n'
#            isotropy_input_str+=' '.join([str(label_index[j]) for j in labels])+'\n'
            isotropy_input_str+=positions
#            print isotropy_input_str
            in_list[centering]=isotropy_input_str
            self.iso_input = isotropy_input_str

        return in_list

    def __get_isotropy_output(self,qe_output=False):
        centering='P'

        in_dict = self.__generate_isotropy_input_from_qe_data(output=qe_output)
        ISODATA = os.path.join(AFLOWpi.__path__[0],'ISOTROPY')
        os.putenv('ISODATA',ISODATA+'/') 

        for centering,in_str in in_dict.iteritems():


            findsym_path = os.path.join(ISODATA,'findsym')
            try:
                find_sym_process = subprocess.Popen(findsym_path,stdin=subprocess.PIPE,stdout=subprocess.PIPE,)

                output = find_sym_process.communicate(input=in_str)[0]
                self.output=output
#                print output



            except:
                print find_sym_process.returncode

            return output

    def __get_sg_num(self):
        sg_info = re.findall('Space Group\s*([0-9]*)\s*([\w-]*)\s*([\w-]*)',self.output)[0]
        self.sgn = int(sg_info[0])

    def __get_iso_cart(self):
            search = 'Lattice vectors in cartesian coordinates:\s*\n(.*\n.*\n.*)\n'
            std_prim_basis_str = re.findall(search,self.output)[0]
            self.iso_pr_car    = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)#*

            return self.iso_pr_car
            


    def __get_iso_basis(self):

        modifier = AFLOWpi.qe.regex.cell_parameters(self.input,return_which='modifier',).lower()
#        if modifier=='angstrom':

#        else:
        self.conv_a = float(re.findall('_cell_length_a\s*([0-9-.]*)',self.output)[0])
        self.conv_b = float(re.findall('_cell_length_b\s*([0-9-.]*)',self.output)[0])
        self.conv_c = float(re.findall('_cell_length_c\s*([0-9-.]*)',self.output)[0])
        self.conv_alpha = float(re.findall('_cell_angle_alpha\s*([0-9-.]*)',self.output)[0])
        self.conv_beta  = float(re.findall('_cell_angle_beta\s*([0-9-.]*)',self.output)[0])
        self.conv_gamma = float(re.findall('_cell_angle_gamma\s*([0-9-.]*)',self.output)[0])
        origin = re.findall('Origin at\s*([0-9-.]*)\s*([0-9-.]*)\s*([0-9-.]*)',self.output)[0]

        origin = numpy.asarray([float(i) for i in origin])
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
                prim_in = AFLOWpi.retr.getCellMatrixFromInput(self.input)
            try:
                prim_in=AFLOWpi.retr._cellStringToMatrix(input_dict['CELL_PARAMETERS']['__content__'])
            except Exception,e:
                print e
                raise SystemExit

        


        self.iso_basis=prim_in



        a=numpy.array([[self.conv_a,],
                       [self.conv_b,],
                       [self.conv_c,],])

#        print self.iso_basis
#        print a
        a=numpy.abs((self.iso_conv).dot(self.iso_basis))

        if self.ibrav in [8,9,10,11]:
            self.conv_a=numpy.sum(a[:,0])
            self.conv_b=numpy.sum(a[:,1])
            self.conv_c=numpy.sum(a[:,2])

#        print self.iso_basis.dot(self.iso_conv)
        
        prim_abc = re.findall('Lattice parameters, a,b,c,alpha,beta,gamma.*\n(.*)\n',self.output)[0].split()
        prim_abc=map(float,prim_abc)
        
        self.axes_flip = numpy.linalg.inv(self.iso_conv)

#        self.axes_flip[0]*=self.conv_a
#        self.axes_flip[1]*=self.conv_b
#        self.axes_flip[2]*=self.conv_c

        self.axes_flip=numpy.around(self.axes_flip,decimals=3)

#        raise SystemExit
        return self.iso_basis

    def __convert_to_ibrav_matrix(self):
        #passing the wrong basis vectors but it's okay since ibrav=ibrav_num overrides it
        #        print self.ibrav


        







        



        
 #       self.qe_basis = AFLOWpi.retr._conv2PrimVec(self.iso_basis,ibrav=self.ibrav)
        self.qe_basis = AFLOWpi.retr.abc2free(a=self.conv_a,b=self.conv_b,c=self.conv_c,alpha=self.conv_alpha,beta=self.conv_beta,gamma=self.conv_gamma,ibrav=self.ibrav,returnString=False)
    


        conv=self.qe_basis*numpy.linalg.inv(self.iso_basis)
#        self.qe_basis[:,0]/=self.conv_a
#        self.qe_basis[:,1]/=self.conv_b
#        self.qe_basis[:,2]/=self.conv_c


        return conv



    def convert(self,ibrav=True):

        input_dict = AFLOWpi.retr._splitInput(self.input)
            
        self.qe_pos = copy.deepcopy(self.orig_pos)


        conv_len = self.axes_flip.T.dot(numpy.array([self.conv_a,self.conv_b,self.conv_c,]).T).tolist()
#        print self.conv_a,self.conv_b,self.conv_c
#        print self.axes_flip
#        self.conv_a,self.conv_b,self.conv_c = numpy.abs(conv_len)

#        self.orig_pos  -= self.origin
        trans = self.iso_basis.dot(numpy.linalg.inv(self.qe_basis))
        self.qe_pos = (self.orig_pos.dot(trans))
        
        



            

#            self.orig_pos[i]-=self.origin
            
#            pos_copy = copy.deepcopy(self.orig_pos[i])
#            pos_copy-=self.origin
#            print pos_copy
#            pre_trans = numpy.matrix(pos_copy)
#            print pre_trans
#            self.qe_pos[i] = self.conv.dot(pos_copy.T).T
#            self.qe_pos[i] = pos_copy.dot(self.conv)
#            pos_copy-=self.origin
#            first = pos_copy.dot(numpy.linalg.inv(self.iso_basis))
            

        qe_pos_str = AFLOWpi.retr._joinMatrixLabels(self.pos_labels,self.qe_pos)
                                         
        modifier = AFLOWpi.qe.regex.cell_parameters(self.input,return_which='modifier',).lower()

        prim_qe_cart = self.iso_pr_car
#        conv_iso = numpy.linalg.inv(self.iso_basis).dot(self.iso_pr_car)

        try:
            del input_dict['CELL_PARAMETERS']
        except:
            pass
            

        input_dict['&system']['a']=self.conv_a#*0.529177249
 #       if self.ibrav in [8,9,10,11,12,13,14]:
        if True:
            input_dict['&system']['b']=self.conv_b#*0.529177249
#        if self.ibrav in [4,6,7,8,9,10,11,12,13,14]:
            input_dict['&system']['c']=self.conv_c#*0.529177249
#        if self.ibrav in [12,13,14]:            
            input_dict['&system']['cosAB']=numpy.cos(self.conv_gamma/180.0*numpy.pi)

#        if self.ibrav in [5]:            
#            input_dict['&system']['cosBC']=numpy.cos(self.conv_gamma/180.0*numpy.pi)
#        if self.ibrav in [14]:                        
            input_dict['&system']['cosAC']=numpy.cos(self.conv_beta/180.0*numpy.pi)
            input_dict['&system']['cosBC']=numpy.cos(self.conv_alpha/180.0*numpy.pi)
#        input_dict['&system']['']=self.conv_a
#        input_dict['&system']['A']=self.conv_a

        input_dict['&system']['ibrav']=self.ibrav

#


        input_dict['ATOMIC_POSITIONS']['__content__']=qe_pos_str

#['__content__']=qe_pos_str
#        input_dict['CELL_PARAMETERS']['__content__']=qe_cm_string
        qe_convention_input = AFLOWpi.retr._joinInput(input_dict)
#        print qe_convention_input
#        print qe_convention_input
        if ibrav==True:
#            print qe_convention_input
            qe_convention_input = AFLOWpi.prep._transformInput(qe_convention_input)
#            print qe_convention_input 
            return qe_convention_input
        else:

            return qe_convention_input




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





    def cif2qe(self):


        input_dict = AFLOWpi.retr._splitInput(self.input)




        convert = AFLOWpi.retr.abc2free(a=1.0,b=1.0,c=1.0,alpha=90.0,beta=90.0,gamma=90.0,ibrav=self.ibrav,returnString=False)

        ins= self.cif.lower()

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

        symops=symops_aux
        ident = numpy.identity(3,dtype=numpy.float64)


        shift,operations = AFLOWpi.prep._process_cif_sym(symops)
        re_atom_pos=re.compile(r'(_atom_site_label.*\n(?:(?:[a-z_\s])*\n))((?:.*\n)*)')
        atom_pos=re_atom_pos.findall(ins)[0]

        loop_list = [x.strip() for x in  atom_pos[0].split('\n') if len(x.strip())!=0]

        spec_lab =  loop_list.index('_atom_site_label')

        x_loc =  loop_list.index('_atom_site_fract_x')
        y_loc =  loop_list.index('_atom_site_fract_y')
        z_loc =  loop_list.index('_atom_site_fract_z')


        positions = [map(str.strip,x.split()) for x in  atom_pos[1].split('\n') if len(x.strip())!=0]

        pos_array=numpy.zeros((len(positions),3))


        labels=[]
        for i in range(len(positions)):
            pos_array[i][0] = float(positions[i][x_loc])
            pos_array[i][1] = float(positions[i][y_loc])
            pos_array[i][2] = float(positions[i][z_loc])

            labels.append(positions[i][spec_lab].strip('0123456789').title())



#



        labels = labels*len(operations)
        all_eq_pos=numpy.zeros((pos_array.shape[0]*operations.shape[0],3))

        for i in xrange(operations.shape[0]):
            temp_pos=numpy.copy(pos_array)
            temp_pos   = temp_pos.dot(operations[i])

            temp_pos[:,0]+=shift[i][0]
            temp_pos[:,1]+=shift[i][1]
            temp_pos[:,2]+=shift[i][2]


            all_eq_pos[i*pos_array.shape[0]:(i+1)*pos_array.shape[0]]=temp_pos

#$        all_eq_pos-=self.origin            
        all_eq_pos%=1.0




        if self.ibrav in [12,13]:
            if numpy.isclose(self.conv_beta,self.conv_gamma):
                all_eq_pos=all_eq_pos[:,[1,2,0]]

            if numpy.isclose(self.conv_alpha,self.conv_gamma):
                all_eq_pos=all_eq_pos[:,[2,0,1]]

            print numpy.linalg.inv(convert.getA())




        b = numpy.ascontiguousarray( all_eq_pos).view(numpy.dtype((numpy.void,  all_eq_pos.dtype.itemsize * all_eq_pos.shape[1])))
        _, idx = numpy.unique(b, return_index=True)

        labels_arr=numpy.array(labels)
        labels_arr=labels_arr[idx]
        all_eq_pos=all_eq_pos[idx]


        if self.ibrav==5:

            convert=numpy.identity(3)

#            beta  = numpy.sqrt((3.0+(self.conv_c/self.conv_a)**2.0))
#            rho_a = beta*self.conv_a/3.0
#            alpha = 2.0*numpy.arcsin(3.0/(2.0*beta))

            convert=self.conv_a*numpy.matrix([[1.0,-1.0*numpy.sqrt(3), 0.0],
                                              [1.0, 1.0*numpy.sqrt(3), 0.0],
                                              [0.0, 0.0           , 2.0*self.conv_c/self.conv_a],])/2.0


            convert=numpy.matrix([[1.0,-1.0*numpy.sqrt(3), 0.0],
                                  [1.0, 1.0*numpy.sqrt(3), 0.0],
                                  [0.0, 0.0           , 2.0],])/2.0



            AH = 1.0#
            CH = 1.0*self.conv_c/self.conv_a
            convert=AH*numpy.matrix([[ 0.0       , AH*numpy.sqrt(3)/1.0, CH],
                                     [ 3.0*AH/2.0,-AH*numpy.sqrt(3)/2.0, CH],
                                     [-3.0*AH/2.0,-AH*numpy.sqrt(3)/2.0, CH],])/3.0
            
#            convert=numpy.linalg.inv(convert)


#        all_eq_pos = numpy.around(all_eq_pos,decimals=10)

        try:
            all_eq_pos=(numpy.linalg.inv(convert.getA()).dot(all_eq_pos.T)).T%1.0
        except:
            all_eq_pos=(numpy.linalg.inv(convert).dot(all_eq_pos.T)).T%1.0

        b = numpy.ascontiguousarray( all_eq_pos).view(numpy.dtype((numpy.void,  all_eq_pos.dtype.itemsize * all_eq_pos.shape[1])))
        _, idx = numpy.unique(b, return_index=True)

        labels_arr=labels_arr[idx]
        all_eq_pos=all_eq_pos[idx]

        spec_sort = numpy.argsort(labels_arr)
        labels_arr=labels_arr[spec_sort]
        all_eq_pos=all_eq_pos[spec_sort]









        atm_pos_str=""
        for i in xrange(all_eq_pos.shape[0]):
            atm_pos_str+= ('%3.3s'% labels_arr[i]) +(' % 9.9f % 9.9f % 9.9f '%tuple(all_eq_pos[i].tolist()))+"\n"



        input_dict['&system']['ibrav']=self.ibrav

        input_dict['&system']['A']=self.conv_a
        input_dict['&system']['B']=self.conv_b
        input_dict['&system']['C']=self.conv_c
        input_dict['&system']['cosAB']=numpy.cos(self.conv_gamma/180.0*numpy.pi)
        input_dict['&system']['cosAC']=numpy.cos(self.conv_beta/180.0*numpy.pi)
        input_dict['&system']['cosBC']=numpy.cos(self.conv_alpha/180.0*numpy.pi)

        if self.ibrav in [4]:
                del input_dict['&system']['cosBC']
                del input_dict['&system']['cosAC']
                del input_dict['&system']['cosAB']
        if self.ibrav in [12,13]:
            if numpy.isclose(self.conv_alpha,self.conv_gamma):
                input_dict['&system']['A']=self.conv_b
                input_dict['&system']['B']=self.conv_c
                input_dict['&system']['C']=self.conv_a
                input_dict['&system']['cosBC']=numpy.cos(self.conv_gamma/180.0*numpy.pi)
                input_dict['&system']['cosAC']=numpy.cos(self.conv_alpha/180.0*numpy.pi)
                input_dict['&system']['cosAB']=numpy.cos(self.conv_beta/180.0*numpy.pi)
            if numpy.isclose(self.conv_beta,self.conv_gamma):
                input_dict['&system']['A']=self.conv_c
                input_dict['&system']['B']=self.conv_a
                input_dict['&system']['C']=self.conv_b
                input_dict['&system']['cosBC']=numpy.cos(self.conv_beta/180.0*numpy.pi)
                input_dict['&system']['cosAC']=numpy.cos(self.conv_gamma/180.0*numpy.pi)
                input_dict['&system']['cosAB']=numpy.cos(self.conv_alpha/180.0*numpy.pi)

    


        if self.ibrav==5:
                    input_dict['&system']['ibrav']=self.ibrav
                    beta  = numpy.sqrt((3.0+(self.conv_c/self.conv_a)**2.0))
                    rho_a = beta*self.conv_a/3.0
                    alpha = 2.0*numpy.arcsin(3.0/(2.0*beta))
                    input_dict['&system']['A']=rho_a
                    input_dict['&system']['B']=rho_a
                    input_dict['&system']['C']=rho_a
                    input_dict['&system']['cosAB']=numpy.cos(alpha)
                    input_dict['&system']['cosAC']=numpy.cos(alpha)
                    input_dict['&system']['cosBC']=numpy.cos(alpha)
        
                
        try:
            del input_dict['CELL_PARAMETERS']
        except:
            pass



        input_dict['ATOMIC_POSITIONS']['__content__']=atm_pos_str
        
        qe_convention_input=AFLOWpi.retr._joinInput(input_dict)
#        if self.ibrav==5:
#            print qe_convention_input
            

        qe_convention_input = AFLOWpi.prep._transformInput(qe_convention_input)
        return qe_convention_input
        #        from scipy.spatial import Delaunay as CH

#         def inside_prim_lat(labels,points,conv,prim):
             
#             labs,ss = AFLOWpi.retr._expandBoundaries(labels,points,2,2,2)
#             ss=ss*2.0-1.0
#  #           prim+=shift
#             prim_hull = numpy.zeros((8,3),dtype=numpy.float64)
# #            prim_hull[0] =
#             prim_hull[1] = prim[0]
#             prim_hull[2] = prim[1]
#             prim_hull[3] = prim[2]
#             prim_hull[4] = prim[0]+prim[1]
#             prim_hull[5] = prim[0]+prim[1]+prim[2]
#             prim_hull[6] = prim[0]+prim[2]
#             prim_hull[7] = prim[1]+prim[2]
            

#             hull=CH(prim_hull)
#             print  hull.find_simplex(ss)
#             in_hull_mask= hull.find_simplex(ss)>=0
#             return labs[in_hull_mask],ss[in_hull_mask]%1.0


#        labels_arr,all_eq_pos= inside_prim_lat(labels_arr,all_eq_pos,ident,convert)
        
