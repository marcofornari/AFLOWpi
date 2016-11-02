import AFLOWpi
import numpy
import os
import subprocess
import StringIO
import re
import copy

class isotropy():

    def __init__(self,input_str,accuracy=0.000,output=False):

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
            

            self.qe   = self.convert()
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
            alat = float(alatSearch.findall(self.input)[-1])
            


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
#            print self.iso_pr_car
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
        iso_conv = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)

#        std_prim_basis = numpy.linalg.inv(AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str))
        std_prim_basis = AFLOWpi.retr._cellStringToMatrix(std_prim_basis_str)
#        print globals().keys()


        input_dict = AFLOWpi.retr._splitInput(self.input)
        if input_dict['CELL_PARAMETERS']["__content__"]=='':
            prim_in = AFLOWpi.retr.getCellMatrixFromInput(self.input)
        else:
            try:
                prim_in=AFLOWpi.retr._cellStringToMatrix(input_dict['CELL_PARAMETERS']['__content__'])
            except Exception,e:
                print e
                raise SystemExit
#$        iso_conv

#        print prim_in*numpy.linalg.inv(iso_conv)

        self.iso_basis=prim_in
#        self.iso_basis=AFLOWpi.retr._cellStringToMatrix(prim_in)

        


#        self.iso_basis[:,0]/=self.conv_a
#        self.iso_basis[:,1]/=self.conv_b
#        self.iso_basis[:,2]/=self.conv_c
        a= numpy.linalg.inv(iso_conv)
        a[:,0]*=self.conv_a
        a[:,1]*=self.conv_b
        a[:,2]*=self.conv_c
        print a
#        self.iso_basis=a
#        self.iso_basis=numpy.around(self.iso_basis,decimals=4)

        self.axes_flip = iso_conv.dot(self.iso_basis)
        self.axes_flip = numpy.linalg.inv(iso_conv)
#        print self.axes_flip
#        self.iso_basis = self.axes_flip.dot(self.iso_basis)
#            self.iso_basis = .dot(self.iso_basis))
#        print numpy.around(numpy.linalg.inv(iso_conv.dot(self.iso_basis)),decimals=5)
#        raise SystemExit

#        self.iso_basis[:,0]/=0.529177249
#        self.iso_basis[:,1]/=0.529177249
#        self.iso_basis[:,2]/=0.529177249
#        print self.iso_basis
#        self.iso_basis-=self.origin
 #       print self.iso_basis
#        raise SystemExit
        return self.iso_basis

    def __convert_to_ibrav_matrix(self):
        #passing the wrong basis vectors but it's okay since ibrav=ibrav_num overrides it
#        print self.ibrav


        

        qe_basis_inv = AFLOWpi.retr._prim2ConvMatrix(self.iso_basis,ibrav=self.ibrav)

        qe_conv2prim = numpy.linalg.inv(qe_basis_inv)
#       self.qe_basis = qe_conv2prim
#        print AFLOWpi.retr.getCellMatrixFromInput(self.input,string=False)

        

#        self.qe_basis =  AFLOWpi.retr.abc2free(self.conv_a,self.conv_b,self.conv_c,self.conv_alpha,self.conv_beta,self.conv_gamma,ibrav=self.ibrav,returnString=False)
#        self.qe_basis =  AFLOWpi.retr.abc2free(self.conv_a,self.conv_b,self.conv_c,self.conv_alpha,self.conv_beta,self.conv_gamma,ibrav=self.ibrav,returnString=False)
        
        self.qe_basis = AFLOWpi.retr._prim2ConvMatrix(self.iso_basis,ibrav=self.ibrav)

        

        self.qe_basis = numpy.linalg.inv(self.qe_basis)
        conv=self.qe_basis*numpy.linalg.inv(self.iso_basis)
        self.qe_basis[:,0]*=self.conv_a
        self.qe_basis[:,1]*=self.conv_b
        self.qe_basis[:,2]*=self.conv_c

        conv = numpy.around(conv,decimals=5)
        return conv
#        prim_in = numpy.round(prim_in/0.529177249,decimals=6)

#        conv = self.qe_basis.dot(numpy.linalg.inv(prim_in))


#        return conv

    def convert(self,ibrav=True):

        input_dict = AFLOWpi.retr._splitInput(self.input)
#        t= numpy.linalg.inv(self.iso_basis)
#        t[:,0] = t[:,0]*self.conv_a*0.529177249
#        t[:,1] = t[:,1]*self.conv_b*0.529177249
#        t[:,2] = t[:,2]*self.conv_c*0.529177249
#        asdf = AFLOWpi.retr._cellStringToMatrix(input_dict['CELL_PARAMETERS']['__content__'])
#        print asdf
#        print t
#        print AFLOWpi.retr._cellMatrixToString(asdf.dot(numpy.linalg.inv(t)))
        
#        print self.qe_basis

#        print self.conv


#        self.conv = conv_conv.dot(self.conv)
#        print self.iso_basis
#        print self.qe_basis

#        print trans
            
        self.qe_pos = copy.deepcopy(self.orig_pos)
#        print self.conv
#        print self.conv
#        self.iso_babsis= self.axes_flip
#        print self.output        
#        self.iso_basis = self.iso_basis.dot()
#        self.conv= self.conv.
        print self.iso_basis
        print self.qe_basis
        conv = self.qe_basis.dot(numpy.linalg.inv(self.iso_basis))
        conv = numpy.around(conv,decimals=4)

        for i in range(len(self.orig_pos)):
            self.origin= numpy.matrix(self.origin)
            mat_pos = numpy.matrix(self.orig_pos[i])
            second = conv.dot(mat_pos.T).T

            self.qe_pos[i]=second.flatten().tolist()[0]
            

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

#        ibrav,cdm1,cdm2,cdm3,cdm4,cdm5,cdm6 =  AFLOWpi.retr.abc2celldm(self.conv_a,self.conv_b,self.conv_c,self.conv_alpha,self.conv_beta,self.conv_gamma,ibrav=self.ibrav,)
#        input_dict['CELL_PARAMETERS']['__content__']= AFLOWpi.retr._cellMatrixToString(
#        input_dict['&system']['celldm(1)']=cdm1
#        input_dict['&system']['celldm(2)']=cdm2
#        input_dict['&system']['celldm(3)']=cdm3
#        input_dict['&system']['celldm(4)']=cdm4
#        input_dict['&system']['celldm(5)']=cdm5
#        input_dict['&system']['celldm(6)']=cdm6


        input_dict['&system']['a']=self.conv_a#*0.529177249
        input_dict['&system']['b']=self.conv_b#*0.529177249
        input_dict['&system']['c']=self.conv_c#*0.529177249
        input_dict['&system']['cosAB']=numpy.cos(self.conv_gamma/180.0*numpy.pi)
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
            qe_convention_input = AFLOWpi.prep._transformInput(qe_convention_input)
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
        elif self.sgn in [143,144,145,146,147,148,149,150,151,152,153,154,
                     155,156,157,158,159,160,161,162,163,164,165,166,
                     167]:
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
#       elif self.sgn in [38,39,40,41]:
#           return 9
#       elif self.sgn in [20,21,35,36,37,63,64,65,66,67,68]:
#           return -9
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

    

#    def __qe_conventional(self):
#        if self.ibrav in [1,2,3,6,7,8,-9,9,10,11]:
#            qe_conv = numpy.asarray([[self.conv_a,0.0,0.0,],
#                                     [0.0,self.conv_b,0.0,],
#                                     [0.0,0.0,self.conv_c,],])

    
#        return qe_conv

