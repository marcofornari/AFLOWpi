import AFLOWpi
import os
import re
import copy 

def __collect_fd_field_forces(oneCalc,ID,for_type='raman'):
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

def __pull_polarization(oneCalc,ID):
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

    e_pol_split = [float(x.split()[-1]) for x in ele_block.split('\n') if len(x.strip())!=0]
    e_ion_split = [float(x.split()[-1]) for x in ion_block.split('\n') if len(x.strip())!=0]
    
    
    e_tot_split = [e_pol_split[0],e_pol_split[1],e_pol_split[2],]

    epol_out_string=''

    epol_out_string+='%18.12f%18.12f%18.12f\n'%(float(e_tot_split[0]),float(e_tot_split[1]),float(e_tot_split[2]),)

    force_postfix='.'.join(ID.split('.')[1:])
    force_postfix=force_postfix.split('_')[0]

    with open(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'epol.%s'%force_postfix),'w+') as out_file_obj:
        out_file_obj.write(epol_out_string)


def __pull_eps_out(oneCalc,ID):
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
    eps_regex=re.compile('Dielectric tensor\n.*\n([0-9\s.-]*)')
    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_epol.out'%ID),'r') as out_file_obj:
            out_string = out_file_obj.read()

        eps_string=eps_regex.findall(out_string)[-1]

        
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    return eps_string

def __pull_born_out(oneCalc,ID):
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
            born_string+='         %s    \n'%num_index
            born_string+=born_charges[i]

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    return born_string

def __gen_fd_input(oneCalc,ID,for_type='raman',de=0.003):
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

    forces = AFLOWpi.run.__collect_fd_field_forces(oneCalc,ID,for_type=for_type)
    header+=forces
    
    finite_field_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_%s.in'%(ID,for_type))
    with open(finite_field_input,'w') as force_file:
        force_file.write(header)

def __swap_scf_inputs(oneCalc,ID):
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

def __setup_raman(calc_subset,orig_oneCalc,field_strength=0.001,nberrycyc=3,for_type='raman'):
    """
    Sets up the input files and the command to run the finite field calculations.
    
    Arguments:
          calc_subset (dict): dictionary of dictionaries of the FD_PHONON calculations 
          orig_oneCalc (dict): oneCalc of one calculation the main calc set

    Keyword Arguments:
         feild_strength (float): applied electric field strength
         nberrycyc (int) number of berry phase cycles
         for_type (str): which type of FD input do we want to generate files for. 'raman','born','eps'

    Returns:
         None
          
    """

    input_string=orig_oneCalc['_AFLOWPI_INPUT_']
    subset_keys = calc_subset.keys()

    in_dict=AFLOWpi.retr.__splitInput(input_string)
    #if possibly a metal don't do LOTO
    if 'degauss' in in_dict['&system'].keys():
        return
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
    
    for index in num:        
        fd_input = AFLOWpi.run.__field_factory(input_string,field_strength=field_strength,for_type=for_type,nberrycyc=nberrycyc,dir_index=index)

        subset_key_index=index%subset_size

        ID = subset_keys[subset_key_index]
        oneCalc=calc_subset[ID]
        
        FD_ID='finite_field.%02d'%index
        file_path_fd_field_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%FD_ID)
    
        print ID

        with open(file_path_fd_field_in,'w') as fd_in_obj:
            fd_in_obj.write(fd_input)


        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN',"""oneCalc = AFLOWpi.run.__swap_scf_inputs(oneCalc,'%s')""" % FD_ID)       
        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN','''AFLOWpi.prep.__to_local_scratch(oneCalc,'%s')'''%FD_ID)

        AFLOWpi.run.__onePrep(oneCalc,ID,engine='espresso',calcType='scf',alt_ID=FD_ID,execPostfix='-northo 1')
        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN',"""AFLOWpi.run.__pull_forces(oneCalc,'%s')""" % FD_ID)       
        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN',"""AFLOWpi.run.__pull_polarization(oneCalc,'%s')""" % FD_ID)       


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

    AFLOWpi.prep.updatelogs(calc_subset_extended,chain_logname,runlocal=True)

def __field_factory(input_str,field_strength=0.001,for_type='raman',nberrycyc=3,dir_index=0):
    """
    Modifies the QE pwscf input for the finite field calc of a given index
    
    Arguments:
         input_str (str): string of QE pwscf input file

    Keyword Arguments:
         feild_strength (float): applied electric field strength
         nberrycyc (int) number of berry phase cycles
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


    in_dict=AFLOWpi.retr.__splitInput(input_str)
    in_dict['&control']['lelfield']='.TRUE.'
    in_dict['&control']['calculation']="'scf'"
    in_dict['&control']['nberrycyc']=nberrycyc
    in_dict['&control']['tprnfor']='.TRUE.'
    in_dict['&control']['prefix']="'_finite_field.%02d'"%dir_index

    input_list=[]

#    for i in range(len(field_dir)):
    for j in range(len(field_dir[dir_index])):
        in_dict['&electrons']['efield_cart(%s)'%(j+1)]=field_dir[dir_index][j]*field_strength

    new_input = AFLOWpi.retr.__joinInput(in_dict)

        
    return new_input

#from __future__ import print_function
import AFLOWpi
import ast
import functools
import os
import subprocess
import glob
import re
import numpy


def __phonon_band_path(oneCalc,ID):
    """
    Gets path for the cell between High Symmetry points in Brillouin Zone
    and returns for to be part of matdyn.x input for the phonon dispersion
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          path (str): Path between High Symmetry points in Brillouin Zone
          
    """

    dk=0.001
    nk=4000
    path = AFLOWpi.retr.__getPath(dk,oneCalc,ID=ID)

    '''scale the k points in the path list so they're as close to nk as possible'''
    splitPath =  [[x.split()[0],int(x.split()[1])] for x in  path.split('\n')[1:] if len(x)!=0]
    total =  [int(x.split()[1]) for x in  path.split('\n')[1:] if len(x)!=0]

    for entries in range(len(total)):
	    if total[entries]==0:
		    total[entries]+=1

    total=sum(total)
    scaleFactor = float(nk)/float(total)

    dk_prime = dk/scaleFactor

    path = AFLOWpi.retr.__getPath(dk_prime,oneCalc,ID=ID)

    return path

def __pull_forces(oneCalc,ID):
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
    

    force_regex=re.compile(r'Forces acting on atoms.*\n\n((?:.*force\s*=\s*[-0-9.]+\s*[-0-9.]+\s*[-0-9.]+\s*)+\n)',re.M)
    try:
        force_block=force_regex.findall(out_string)[-1]
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
    force_split = [x.split()[6:] for x in force_block.split('\n') if len(x.strip())!=0]

    force_out_string=''
    for i in force_split:
        force_out_string+='%18.12f%18.12f%18.12f\n'%(float(i[0]),float(i[1]),float(i[2]),)

    force_postfix='.'.join(ID.split('.')[1:])
    force_postfix=force_postfix.split('_')[0]

    with open(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'force.%s'%force_postfix),'w+') as out_file_obj:
        out_string = out_file_obj.write(force_out_string)



import os
def __pp_phonon(__submitNodeName__,oneCalc,ID,LOTO=True,de=0.01,raman=True,field_strength=0.001):
    """
    Calls the executables for the post processing of the force
    data generated by the scf calculations created by fd.x
    
    Arguments:
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          None
          
    """
        
    if AFLOWpi.prep.__ConfigSectionMap("run","execprefix") != '':
	    execPrefix=AFLOWpi.prep.__ConfigSectionMap("run","execprefix")
    else:
        execPrefix=''
            
    #run fd_ifc
    AFLOWpi.run.__oneRun(__submitNodeName__,oneCalc,'%s_fd_ifc'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./fd_ifc.x' ) 

    # get ifc file for matdyn
    header_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'FD_PHONON','header.txt')
    with open(header_file,'r') as header_fileobj:
        header_string=header_fileobj.read()
        


    if raman==True:
        field_list=['epol','born','raman']
    elif LOTO==True:
        field_list=['epol','born',]
    else:
        field_list=[]

    try:
        for field in field_list:
            try:
                AFLOWpi.run.__gen_fd_input(oneCalc,ID,for_type=field,de=field_strength)
                AFLOWpi.run.__oneRun(__submitNodeName__,oneCalc,'%s_%s'%(ID,field),execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./fd_ef.x' )


            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)            

#######################################################################
        if raman==True or LOTO==True:
            epol_string=AFLOWpi.run.__pull_eps_out(oneCalc,ID)
            born_charge_string=AFLOWpi.run.__pull_born_out(oneCalc,ID)

            LOTO_REPLACE='''T

    %s
    %s

    ''' % (epol_string,born_charge_string)


            header_string=re.sub("F\n",LOTO_REPLACE,header_string)
########################################################################
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        
    
    ifc_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'PHONON_ifc.fc')
    with open(ifc_file,'r') as ifc_fileobj:
        ifc_string=ifc_fileobj.read()

    for_matdyn  = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.fc' % ID)
    with open(for_matdyn,'w') as mat_dyn_in_fileobj:
        mat_dyn_in_fileobj.write(header_string)
        mat_dyn_in_fileobj.write(ifc_string)


    #run matdyn for bands
    AFLOWpi.run.__oneRun(__submitNodeName__,oneCalc,'%s_matdyn_phBand'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )

    #run matdyn for dos
    AFLOWpi.run.__oneRun(__submitNodeName__,oneCalc,'%s_matdyn_phDOS'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )
    
    #remove the fd.x output files in case we do another phonon calc
    try:
	    globpath=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'FD_PHONON')
	    phil = glob.glob(globpath+'/displaced*.in')
	    for i in phil:
		    os.remove(i)
    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)


def write_fdx_template(oneCalc,ID,nrx1=2,nrx2=2,nrx3=2,innx=2,de=0.01,atom_sym=True,disp_sym=True):
    """
    Generates input files for fd.x, fd_ifc.x, and matdyn.x for the finite
    difference phonon calculations. 
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          nrx1 (int): supercell size for first primitive lattice vector
	  nrx2 (int): supercell size for second primitive lattice vector
	  nrx3 (int): supercell size for third primitive lattice vector
	  innx (int): how many differernt shifts in each direction for 
                      finite difference phonon calculation
	  de (float): amount to shift the atoms for finite differences

    Returns:
          None          
          
    """


    if atom_sym==True:
        noatsym='.false.'
    else:
        noatsym='.true.'
    if disp_sym==True:
        nodispsym='.false.'
    else:
        nodispsym='.true.'

    temp_dir= AFLOWpi.prep.__get_tempdir()

    fd_template='''!fd_prefix     Prefix of the preceding pw.x run
!fd_outdir     Outdir of the preceding pw.x run
!fd_outfile    Prefix for the generated macrocell input files
!fd_outdir_dir Directory for the generated macrocell input files
!nrx1,2,3      Number of unitcell repetitions along the lattice vectors
!              (for the generation of the macrocell)
!de            Cartesian displacement for the atoms, in Angs.
!
&inputfd
 fd_prefix      = '%s' 
 fd_outdir      = './'
 fd_outfile     = 'displaced'
 fd_outfile_dir = './FD_PHONON'
 noatsym        = %s
 nodispsym      = %s
 nrx1           = %i
 nrx2           = %i
 nrx3           = %i
 innx           = %i
 de             = %f
/

!The following namelist is for additional flexibility
!1. Each string (everything between a pair of single quotation marks)
!   will be pasted verbatim to the generated input files.
!   Within the corresponding namelist
!2. Mind not using single quotation marks within strings
!3. Mind keeping "system2" instead of "system". Fortran conflict
!   within the corresponding
!4. Asterisks (*) inside the "kpoints" string represent change of line
!
&verbatim
 control =
'
'

 electrons =
'
'

 system2 =
'
'

 kpoints =
'
'
/
'''%(oneCalc['_AFLOWPI_PREFIX_'],noatsym,nodispsym,nrx1,nrx2,nrx3,innx,de)
    FD_in_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_fd.in'%ID)
    with open(FD_in_path,'w') as fd_in_file:
        fd_in_file.write(fd_template)



    fd_ifc_template='''&input
      prefix='%s'
!      outdir='./'
      nrx1=%i,
      nrx2=%i,
      nrx3=%i,
      innx=%i,
      de=%f,
      file_force='./FD_PHONON/force',
      file_out='./PHONON_ifc'
      noatsym        = %s
      nodispsym      = %s
      verbose=.true.
      hex=.false.
    /
     '''%(oneCalc['_AFLOWPI_PREFIX_'],nrx1,nrx2,nrx3,innx,de,noatsym,nodispsym)

    FD_ifc_in_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_fd_ifc.in'%ID)
    with open(FD_ifc_in_path,'w') as fd_ifc_in_file:
        fd_ifc_in_file.write(fd_ifc_template)




    inputDict=AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])

    masses=[]
    species_string = inputDict['ATOMIC_SPECIES']['__content__']
    
    split_species_string = species_string.split('\n')
    for i in split_species_string:
        try:
            masses.append(float(i.split()[1]))
        except:
            pass

    amass_str=''
    for i in range(len(masses)):
        amass_str+='   amass(%d) = %f,\n'% (i+1,masses[i])


    ph_band_path=__phonon_band_path(oneCalc,ID)



    matdyn_template='''
 &input
   asr='crystal',
%s
   flvec='%s.modes'
   flfrc='%s.fc',
   flfrq='%s.phBAND', 
   fldyn='%s.dyn',
   fleig='%s.eig',
!   l1=%d
!   l2=%d
!   l3=%d
   nosym=%s
   fd=.true.
   na_ifc=.true.
   q_in_band_form=.true.,
 /
%s
'''%(amass_str,ID,ID,ID,ID,ID,nrx1,nrx2,nrx3,noatsym,ph_band_path)
    matdyn_in_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_matdyn_phBand.in'%ID)
    with open(matdyn_in_path,'w') as matdyn_in_file:
        matdyn_in_file.write(matdyn_template)


    kgrid=''
    kpt_str  = inputDict['K_POINTS']['__content__']    
    kpt_ints = [int(numpy.ceil(float(x)*5)) for x in kpt_str.split()[:3]]

    for i in range(len(kpt_ints)):
        kgrid+='   nk%d = %s,\n'%(i+1,kpt_ints[i])

    matdyn_dos_template='''
 &input
   asr='crystal',
%s
   fd=.true.
   na_ifc=.true.
   fldos='%s.phdos'
   flfrc='%s.fc',
   flfrq='%s.phDOS', 
   nosym=%s
!   l1=%d
!   l2=%d
!   l3=%d
   dos=.true.,

%s   
 /

'''%(amass_str,ID,ID,ID,noatsym,nrx1,nrx2,nrx3,kgrid)
    matdyn_in_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_matdyn_phDOS.in'%ID)
    with open(matdyn_in_path,'w') as matdyn_in_file:
        matdyn_in_file.write(matdyn_dos_template)



def clean_cell_params(output):
    """
    Parses the atomic shifts in a supercell from fd.x
    outputted pw.x input files to correct for formatting
    issues when they are imported to AFLOWpi
    
    Arguments:
          output (str): pw.x input files generated by fd.x

    Keyword Arguments:
          None

    Returns:
          output (str): pw.x input files generated by fd.x (cleaned by AFLOWpi)
          
    """

    lines = iter(output.splitlines())
    output=''
    for line in lines:
        # signifies start of shifts in supercell
        if 'CELL_PARAMETERS' not in line:
            output+=line+'\n'
            continue
        else:
            break
        # use same iterator so we don't have to store an index
    params = [[ast.literal_eval(n) for n in line.split()] for line in lines]
    output+='CELL_PARAMETERS {angstrom}\n'
    for i in range(len(params)):    
	    for j in range(len(params[i])):
		    if numpy.abs(params[i][j])<0.0001:
			    params[i][j]=0.0
			    

    for j in params:
        output+=' '.join([str(i) for i in j])+'\n'


    return output


def prep_fd(__submitNodeName__,oneCalc,ID,nrx1=2,nrx2=2,nrx3=2,innx=2,de=0.01,atom_sym=True,disp_sym=True):

    """
    Generates input files for fd.x, fd_ifc.x, and matdyn.x for the finite
    difference phonon calculations. 
    
    Arguments:
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          nrx1 (int): supercell size for first primitive lattice vector
	  nrx2 (int): supercell size for second primitive lattice vector
	  nrx3 (int): supercell size for third primitive lattice vector
	  innx (int): how many differernt shifts in each direction for 
                      finite difference phonon calculation
	  de (float): amount to shift the atoms for finite differences

    Returns:
          None          
          
    """

    #copy fd.x to the directories
    engineDir  = AFLOWpi.prep.__ConfigSectionMap("prep",'enginedir')

    for fd_exec in ['fd.x','fd_ifc.x','fd_ef.x','matdyn.x']:
        if AFLOWpi.prep.__ConfigSectionMap('prep','copyexecs').lower()!='false':
            AFLOWpi.prep.totree(os.path.join(engineDir,'%s'%fd_exec),{ID:oneCalc},symlink=False)
        else:
            AFLOWpi.prep.totree(os.path.join(engineDir,'%s'%fd_exec),{ID:oneCalc},symlink=True)

    write_fdx_template(oneCalc,ID,nrx1=nrx1,nrx2=nrx2,nrx3=nrx3,innx=innx,de=de,atom_sym=atom_sym,disp_sym=disp_sym)
    


    AFLOWpi.run.__oneRun(__submitNodeName__,oneCalc,'%s_fd'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./fd.x' )    


    ocd = oneCalc['_AFLOWPI_FOLDER_']
    globpath=os.path.join(ocd,'FD_PHONON')

    phil = glob.glob(globpath+'/displaced*.in')


    for i in phil:
        fdfo = open(i,'r')
        fdfos = fdfo.read()
        fdfo.close()
        NC=clean_cell_params(fdfos)
        fdfo = open(i,'w')
        fdfo.write(NC)
        fdfo.close()
        
#    return fd_pwscf_in
    # pass repr of input_ to ensure it's fully quoted
#    cmd = 'fd.x < {0!r}'.format(input_)

#    run = monitor('./fd_files')(subprocess.check_output)

#    changed, output = run(cmd, shell=True, universal_newlines=True)

#    return changed, output



# if __name__ == '__main__':
#     paths, output = fd('fd.in')

#     for path in paths:
        
#         with open(path) as f:
            
#             params = cell_params(f.read())

# if params:
#     print(path, params)



#    return fd_template
# def monitor(dirpath):
#     """Decorator that monitors `dirpath` non-recursively and returns a set of files changed by the decorated function."""

#     # avoid issues with relative paths
#     abspath = os.path.abspath(dirpath)

#     # mappable callable to again avoid issues with relative paths for changed files
#     joiner = lambda f: os.path.join(abspath, f)

#     def wrapper(func):
#         @functools.wraps(func)

# def wrapped(*args, **kwargs):
#     try:
#         # time of last access
#         before = {f: os.path.getmtime(f) for f in map(joiner, os.listdir(abspath))}

#     # dirpath does not exist yet, any files added must be from function call
#     except os.error:
#         before = {}

#     result = func(*args, **kwargs)

#     # use sets to avoid implying that files were changed in a specific order
#     try:
#         # new files won't be in before dict
#         changed = {f for f in map(joiner, os.listdir(abspath)) if not (os.path.getmtime(f) == before.get(f))}

#     # dirpath still does not exist, best not to raise an error
#     except os.error:
#         changed = set()
    
#     return changed, result

# return wrapped

# return wrapper

import os
import datetime
import cPickle
import logging 
import re 
import subprocess
#import AFLOWpi.prep
import sys                                                                     
import signal
import collections
import copy
import logging.handlers
import AFLOWpi
import time
import Queue
import socket
import inspect
import traceback
import random
import string
import atexit
import glob
import __main__
import numpy
import AFLOWpi
import itertools

if __name__!='__main__':
    engineDict={}

def __exitClean(signal,frame):
    """
    If the terminate 15 signal is received exit with 0 exit code
    
    Arguments:
          signal (signal.SIGNAL): *nix signal
          frame (frame): Inspect frame

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    logging.debug('got the signal to exit cleanly. Exiting.')
    sys.exit(0)

def __recordDeath(signal,frame):
    """
    If the 15 signal to terminate the process is sent. 
    Record it in oneCalc['__status__']['Error']
    
    Arguments:
          signal (signal.SIGNAL): *nix signal
          frame (frame): Inspect frame

    Keyword Arguments:
          None

    Returns:
          None          
          
    """

    logging.debug('Job Killed.')

    index=globals()["__TEMP__INDEX__COUNTER__"]
    if int(__main__.oneCalc['_AFLOWPI_INDEX_'])==index:
        __main__.oneCalc['__status__']['Error']="Killed"
        AFLOWpi.prep.__saveOneCalc(__main__.oneCalc,__main__.ID)
    
def cleanup(calcs):
    """
    Deletes all files a calculation set's  directory
    tree that are prepended with '_'
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations                

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.prep.__addToBlock(oneCalc,ID,'CLEANUP','\nimport os\nos.system("rm -rf %s/_*")\n' % oneCalc['_AFLOWPI_FOLDER_'])            
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)



def addatexit__(command,*args,**kwargs):
    """
    Wrapper to add function to be run at exit
    
    Arguments:
          command (func): function to be run
          *args: arguments for command

    Keyword Arguments:
          **kwargs: Keyword arguments for command

    Returns:
          None
          
    """
    
    if '__atexitList' not in globals().keys():
        #######################################
        @atexit.register
        def __runatexitList():
            for item in __atexitList:
                item[0](*item[1],**item[2])
        global __atexitList
        __atexitList = []

    __atexitList.append((command,args,kwargs))




def __colorize_message(string,level='ERROR',show_level=True):
    """
    Colorizes text. Used for colorizing the logging
    
    Arguments:
          string (str): A string of text

    Keyword Arguments:
          level (str): Specific colors are chosen for logging message type

    Returns:
          levelname_color (str): Colorized version of the string input
          
    """

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[1;%dm"
    BOLD_SEQ = "\033[1m"

    COLORS = {
        'WARNING': YELLOW,
        'INFO': WHITE,
        'DEBUG': BLUE,
        'CRITICAL': YELLOW,
        'ERROR': RED
        
    }



    
    levelname=level
    if show_level==True:
        levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname +": "+string+ RESET_SEQ
    else:
        levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + string+ RESET_SEQ
    return levelname_color







def _fancy_error_log(e):
    """
    Logs an error and prints it out on the stdout if logging=debug in the config file used
    
    Arguments:
          e (str): string of the error 

    Keyword Arguments:
          None

    Returns:
          None
          
    """

#    e = __colorize_message(e)

    logging.error(__colorize_message('%%%DEBUG%%%'))
    logging.error(e)
    _, _, tb = sys.exc_info()
    errorList =  traceback.format_list(traceback.extract_tb(tb)[-6:])[-6:]
    for errorMSG in errorList:
        logging.error(__colorize_message(errorMSG))     
        logging.error(__colorize_message('%%%DEBUG%%%'))       

        
    if AFLOWpi.prep.__ConfigSectionMap('prep','loglevel').upper() == 'DEBUG':
        print '%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%'
        try:
            print e
            for errorMSG in errorList:
                print errorMSG

        except Exception,e:
            print e

        print '%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%'



    


def  __getEnginePath(engine,calcType):
    """
    Gives the name of the executable file for the ab initio engine
    for a given type of calculation.
    
    Arguments:
          engine (str): Ab initio engine being used
          calcType (str): Type of calculation to be done

    Keyword Arguments:
          None

    Returns:
          execPath (str): the name of the executable for that engine for that calcType
          
    """

    if engine=='espresso':
        execDict={'bands':['./bands.x'],

                  'pdos':['./projwfc.x'],
                  'dos':['./dos.x'],
                  'scf':['./pw.x'],}
    else:
        execDict={}
    try:
        execPath=execDict[calcType]
    except:
#        print 'Can not find executable path. May not be implemented Returning Blank String'
 #       logging.debug('Can not find executable path. May not be implemented. Returning Blank String')
        execPath=''
    return execPath

################################################################################################################

################################################################################################################
def  __getExecutable(engine,calcType):
    """
    Gives the name of the executable file for the ab initio engine
    for a given type of calculation.

    OBSOLETE. NEEDS REMOVAL. AFLOWpi.run.__getEnginePath is almost identical 
    
    Arguments:
          engine (str): Ab initio engine being used
          calcType (str): Type of calculation to be done

    Keyword Arguments:
          None

    Returns:
          executable (str): the name of the executable for that engine for that calcType
          
    """



    if engine=='espresso':
        execDict={'scf':'pw.x',
                  'dos':'dos.x',
                  'pdos':'projwfc.x',
                  'bands':'bands.x',
                  'emr':'gipaw.x', 
                  'nmr':'gipaw.x', 
                  'hyperfine':'gipaw.x', 
                  'gvectors':'gipaw.x',
                  }

    elif engine.lower()=='want':
        execDict={'bands','bands.x',}

    else:
        execDict={}
    try:
        executable=execDict[calcType]
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
#        print 'Can not find executable path. May not be implemented Returning Blank String'
#        logging.debug('Can not find executable path. May not be implemented. Returning Blank String')
        execPath=''
    return executable




import tempfile
import pipes
################################################################################################################

################################################################################################################

# def __fixconverge(ID,folder):
#     logging.debug('entering __fixconverge')
#     logging.warning('problem converging %s/%s.in' % (folder,ID))
#     logging.debug('exiting __fixconverge') 


################################################################################################################


# def __resubmit(ID,oneCalc,__submitNodeName__,queue,calcName):

#         def tryAgain():
#             __resubmit(ID,oneCalc,__submitNodeName__,queue,calcName)

#         signal.signal(signal.SIGTERM, tryAgain)	
#         try:

#             IDsplit = ID.split('_')[0]
#         except:
#             IDsplit=ID

#         newqsubFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.qsub' % ID)

#         print 'resubmitting %s' % calcName
#         command = "ssh -o StrictHostKeyChecking=no %s 'qsub %s -N %s %s'" % (__submitNodeName__,queue,calcName,newqsubFileName)

#         try:
#             subprocess.check_call(command,stdout=subprocess.PIPE,shell=True)            
#             sys.exit(0)


#             print 'resubmission of %s successful' % calcName

#         except Exception,e:
#             AFLOWpi.run._fancy_error_log(e)
#             __resubmit(IDsplit,oneCalc,__submitNodeName__,queue,calcName)


################################################################################################################

def generateSubRef(qsubRefFileString, oneCalc,ID):
    """
    Reads in the reference cluster submission file specified in "jobreffile"
    in the config used. Tries to insert a few parameters.

    OBSOLETE PLANNED FOR REMOVAL

    Arguments:
          qsubRefFileString (str): string of the "reference" cluster submission file
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          clusterTypeDict (dict): string cluster submission file
          
    """



    logging.debug('entering generateSubRef')
    clusterTypeDict={}

    calcName=AFLOWpi.run.__get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])

    stderrFileNm=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_Cluster_'+ID+'.stderr')
    stdoutFileNm=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_Cluster_'+ID+'.stdout')

    qsubRef = re.sub('_AFLOWPI_CALCNM_',calcName,qsubRefFileString)
    qsubRef = re.sub('_AFLOWPI_STDOUT_',stdoutFileNm,qsubRef)
    qsubRef = re.sub('_AFLOWPI_STDERR_',stderrFileNm,qsubRef)

    clusterTypeDict['qsubRef']=qsubRef
    logging.debug('exiting generateSubRef')


    return clusterTypeDict



#############################################################################################################

#########################################################################################################################
def emr(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
    """
    Wrapper to set up GIPAW EMR calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                            (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL
          

    Returns:
          None 
          
    """

    testOne(calcs,calcType='emr',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    gipawdir=AFLOWpi.prep.__ConfigSectionMap('prep','gipawdir')
    if AFLOWpi.prep.__ConfigSectionMap('prep','copyexecs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)

    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='emr')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='emr')
        except Exception,e:
            _fancy_error_log(e)


def hyperfine(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None,isotope=()):
    """
    Wrapper to set up GIPAW hyperfine calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                             (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL
          

    Returns:
          None
          
    """

    testOne(calcs,calcType='hyperfine',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)

    gipawdir=AFLOWpi.prep.__ConfigSectionMap('prep','gipawdir')
    if AFLOWpi.prep.__ConfigSectionMap('prep','copyexecs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)


    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='hyperfine')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='hyperfine')
        except Exception,e:
            _fancy_error_log(e)



def nmr(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
    """
    Wrapper to set up GIPAW NMR calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations               

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                             (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL
          

    Returns:
          None
          
    """

    testOne(calcs,calcType='nmr',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    gipawdir=AFLOWpi.prep.__ConfigSectionMap('prep','gipawdir')
    if AFLOWpi.prep.__ConfigSectionMap('prep','copyexecs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)



    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='nmr')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='nmr')
        except Exception,e:
            _fancy_error_log(e)




def gvectors(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
    """
    Wrapper to set up GIPAW gvectors calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                             (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL


    Returns:
          None          
          
    """

    testOne(calcs,calcType='gvectors',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)

    gipawdir=AFLOWpi.prep.__ConfigSectionMap('prep','gipawdir')
    if AFLOWpi.prep.__ConfigSectionMap('prep','copyexecs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)


    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='gvectors')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='gvectors')
        except Exception,e:
            _fancy_error_log(e)




#########################################################################################################################



#############################################################################################################
def __skeletonRun(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
    """
    Wrapper to set up a custom calculation. Inputs and oneCalc calculation
    dictionary must be created before calling this.
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                             (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL
          
    Returns:
          None
          
    """

    testOne(calcs,calcType=None,engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)

def scf(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
    """
    Wrapper to set up self-consitent calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations              

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                             (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL

    Returns:
          
          
    """

    testOne(calcs,calcType='scf',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf')
        except Exception,e:
            _fancy_error_log(e)


def dos(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
    """
    Wrapper to set up DOS nscf calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                            (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL    

    Returns:
          None
          
    """

    testOne(calcs,calcType='dos',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='dos')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='dos')
        except Exception,e:
            _fancy_error_log(e)

def pdos(calcs,engine='',execPrefix=None,execPostfix='',holdFlag=True,config=None):
    """
    Wrapper to set up DOS projection calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations                    

    Keyword Arguments:
          engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                            (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL

    Returns:
          None
          
    """    

    testOne(calcs,calcType='pdos',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='pdos')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='pdos')
        except Exception,e:
            _fancy_error_log(e)



def bands(calcs,engine='',execPrefix=None,execPostfix=' ',holdFlag=True,config=None):
    """
    Wrapper to set up Electronic Band Structure calculation
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations
          execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)          
          execPostfix (str): commands to go after the executable when run 
                            (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          holdFlag (bool): DEFUNCT. NEEDS REMOVAL
          config (str): DEFUNCT. NEEDS REMOVAL
          

    Returns:
          None
          
    """

    if engine=='':
        engine = AFLOWpi.prep.__ConfigSectionMap("run",'engine')
    
    testOne(calcs,calcType='bands',engine=engine,execPrefix=execPrefix,execPostfix='',holdFlag=holdFlag,config=config)
    for ID,oneCalc in calcs.iteritems():
        try:
            __onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='bands')
        except Exception,e:
            _fancy_error_log(e)
        try:
            __testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='bands')
        except Exception,e:
            _fancy_error_log(e)

        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN','AFLOWpi.run.__PW_bands_fix(oneCalc,ID)')


#############################################################################################################

#############################################################################################################

def testOne(calcs,calcType='scf',engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
        """
        Run all the calculation in the dictionary with a specific engine
    
        Arguments:
               calcs (dict): Dictionary of dictionaries of calculations

        Keyword Arguments:
	       engine (str): executable that you are calling to run the calculations
               execPrefix (str): commands to go before the executable when run 
                                 (ex. mpiexec nice -n 19 <executable>) (default = None)          
               execPostfix (str): commands to go after the executable when run 
                                  (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)

        Returns:
               None          
          
        """
        #############################################################################################################

        try:
            logging.debug('entering testOne')


                
            if execPrefix == None:
                if AFLOWpi.prep.__ConfigSectionMap("run","execprefix") != '':
                    execPrefix=AFLOWpi.prep.__ConfigSectionMap("run","execprefix")

                else:
                    execPrefix=''

            if execPostfix == None:
                if AFLOWpi.prep.__ConfigSectionMap("run","execpostfix") != '':
                    execPostfix = AFLOWpi.prep.__ConfigSectionMap("run","execpostfix")

                else:
                    execPostfix=''

            
            if calcType=='bands':
                if len(re.findall(r'pool',execPostfix))!=0:
                    execPostfix=' '

            if engine=='':
                engine = AFLOWpi.prep.__ConfigSectionMap("run",'engine')

            engineDir  = AFLOWpi.prep.__ConfigSectionMap("prep",'enginedir')	
            if os.path.isabs(engineDir) == False:
                configFileLocation = AFLOWpi.prep.__getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                engineDir =  os.path.join(configFileLocation, engineDir)



            try:
                enginePath =  AFLOWpi.run.__getEnginePath(engine,calcType)
            except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)

            for files in enginePath:

                if os.path.exists(os.path.join(engineDir,files))==False:
                    print '%s not found. Check your config file to make sure enginedir path is correct and that %s is in that directory..Exiting' %(files,files)
                    logging.error('%s not found. Check your config file to make sure enginedir path is correct and that %s is in that directory..Exiting' %(files,files))
                    raise SystemExit

                if AFLOWpi.prep.__ConfigSectionMap('prep','copyexecs').lower()!='false':
                    AFLOWpi.prep.totree(os.path.abspath(os.path.join(engineDir,files)),calcs)
                else:
                    for ID,oneCalc in calcs.iteritems():
                        try:
                            os.symlink(os.path.abspath(os.path.join(engineDir,files)),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],files))
                        except OSError:
                            os.unlink(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],files))
                            os.symlink(os.path.abspath(os.path.join(engineDir,files)),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],files))

            logfile = os.path.abspath(os.path.join((os.path.dirname(calcs[random.choice(calcs.keys())]['_AFLOWPI_FOLDER_'])),'AFLOWpi','LOG.log'))
        except Exception,e:
            _fancy_error_log(e)


        
        #############################################################################################################
        try:
            calcsCopy=copy.deepcopy(calcs)
            newCalcs = collections.OrderedDict()

            for ID,oneCalc in calcs.iteritems():

                calcs[ID]['_AFLOWPI_CONFIG_']=AFLOWpi.prep.__getConfigFile()
                calcs[ID]['calcType']=calcType                

        except Exception,e:
            _fancy_error_log(e)
        #############################################################################################################
        clusterType = AFLOWpi.prep.__ConfigSectionMap("cluster",'type')
        if clusterType=='':
            clusterType=''
        try:
            if 'firstCalcList' not in globals().keys() or holdFlag==False:
#                global firstCalcList
#                firstCalcList = True
                try:
                    baseID = ID.split('_')[0]
                    try:
                        #splitIndexInt = int(splitIndex) 
                        splitIndexInt = int(oneCalc['__TEMP__INDEX__COUNTER__'])
                        prev_ID = '%s_%.02d' % (baseID.split('_')[0],splitIndexInt-1)
                    except:
                        prev_ID=ID
                        pass                   

                    if os.path.exists(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+prev_ID+'.py')):
                        pass

                    else:
                        prev_ID=ID

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)

#                __addatexit(submitFirstCalcs__,calcs)
            elif firstCalcList!=False:
                pass


            firstCalcList=False
        except Exception,e:
            _fancy_error_log(e)

        logging.debug('exiting testOne')


        #############################################################################################################
def reset_logs(calcs):
    """
    Removes log files from AFLOWpi directory
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    for ID,oneCalc in calcs.iteritems():
        if oneCalc['__chain_index__']!=1:
            log_file=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.oneCalc' % ID)            
            os.remove(log_file)
        else:
            oneCalc['__execCounter__']=0
            AFLOWpi.prep.__saveOneCalc(oneCalc,ID)



def resubmit(calcs):
    """
    Stages loaded calculation set to be resubmitted on a cluster
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations

    Keyword Arguments:
          None    

    Returns:
          None
          
    """

    AFLOWpi.run.addatexit__(AFLOWpi.run.submitFirstCalcs__,calcs)
    AFLOWpi.run.submit()


def __get_qsub_name(path):
    """
    Takes path of the cluster submission file and forms a name for the submission
    
    Arguments:
          path (str): path of the cluster job submission file
          
    Keyword Arguments:
          None

    Returns:
          calcName (str): name of the calculation used with the submission
          
    """

    c_name = os.path.basename(os.path.normpath(path))
#    if c_name[-2].strip=='':
#        try:
#            calc_name='_'.join([c_name[-3],c_name[-1]])
#        except:
#            calcName = '_'.join(c_name[-2:])
#    else:
    try:
        calcName = '_'.join([j for j in c_name.split('_')[-4:] if j!=''])
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        try:
            calcName = '_'.join([j for j in c_name.split('_')[-3:] if j!=''])
        except:
            calcName = '_'.join(c_name.split('_')[-2:])

    try:
        calcName=calcName[-16:]
    except:
        pass 

    return calcName



import datetime
def __qsubGen(oneCalc,ID):
    """
    Generates the cluster job submission file
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          qsubFileString (str): string of cluster submission file
          
    """

    logging.debug('ENTERING __qsubGen')

    calcName=AFLOWpi.run.__get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])

    execFileString = os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py'))

    tmpdir_envar=''

    ls_option = AFLOWpi.prep.__ConfigSectionMap("cluster",'localscratch').strip().lower()
    if ls_option=='true':
        tmpdir = AFLOWpi.prep.__ConfigSectionMap("cluster",'localscratchdir').strip()
        if tmpdir!='':
            tmpdir_envar="export TMPDIR='%s'\n"%tmpdir


    try:
        qsubFilename = os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.qsub'))
        prevQsubString = ''
        clusterType = AFLOWpi.prep.__ConfigSectionMap("cluster",'type').strip().upper()
        if clusterType=='PBS':
            qsubRefFileName = AFLOWpi.prep.__ConfigSectionMap("cluster",'jobreffile')
            if os.path.isabs(qsubRefFileName) == False:
                configFileLocation = AFLOWpi.prep.__getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                qsubRefFileName =  os.path.join(configFileLocation, qsubRefFileName)        
            qsubRefString = file(qsubRefFileName,'r').read()
            qsubRef = generateSubRef(qsubRefString,oneCalc,ID)['qsubRef']
            qsubRef = AFLOWpi.prep.remove_blank_lines(qsubRef)

            with open(qsubFilename,'w') as qsubFile:
                qsubSub='''%scd %s
python %s''' % (tmpdir_envar,os.path.abspath(oneCalc['_AFLOWPI_FOLDER_']),os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py')))


                #redirect the -e, -o or -oe to the directory for that job
                dash_e_regex=re.compile(r'#PBS\s+-e\s+.*\n')
                dash_o_regex=re.compile(r'#PBS\s+-o\s+.*\n')
                dash_oe_regex=re.compile(r'#PBS\s+-oe\s+.*\n')

                stderr_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_Cluster_'+ID+'.stderr')
                stdout_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_Cluster_'+ID+'.stdout')
                if len(dash_e_regex.findall(qsubRef)):
                    qsubRef=dash_e_regex.sub('\n#PBS -e %s\n'%stderr_name,qsubRef)
                else:
                    qsubRef+='\n#PBS -e %s\n'%stderr_name
                if len(dash_o_regex.findall(qsubRef)):
                    qsubRef=dash_o_regex.sub('\n#PBS -o %s\n'%stdout_name,qsubRef)
                else:
                    qsubRef+='\n#PBS -o %s\n'%stdout_name
                if len(dash_oe_regex.findall(qsubRef)):
                    qsubRef=dash_oe_regex.sub('\n',qsubRef)
                else:
                    pass

                #put the name of the AFLOWpi job in the  PBS submit script

                name_regex=re.compile(r'\s*#PBS.*-[nN].*\n')
                calcNameString = '\n#PBS -N %s\n'%calcName
                if len(name_regex.findall(qsubRef)):
                    qsubRef=name_regex.sub(calcNameString,qsubRef)
                else:
                    qsubRef=calcNameString+qsubRef
                if len(re.findall('_AFLOWPI_QSUB_',qsubRef)):
                    qsubRef = re.sub('_AFLOWPI_QSUB_',qsubSub,qsubRef)
                else:
                    qsubRef+='\n'+qsubSub+'\n'


                qsubFileString=qsubRef
                qsubFile.write(qsubFileString)

    except Exception,e:
        _fancy_error_log(e)
    try:
            '''
            record the walltime requested and whether or not user is using stepsasjobs=false
            in a file so that if we need to do a restart we can figure out how much time we 
            have left to run for the initial submission with stepsasjobs=false
            '''
#            walltime=__getWalltime(oneCalc,ID)
#            stepsasjobs = AFLOWpi.prep.__ConfigSectionMap("cluster",'stepsasjobs').strip().upper()
#            if stepsasjobs.lower() != 'false':
#                stepsasjobs='true'
            pass
    except Exception,e:
        _fancy_error_log(e)

    return qsubFileString


#############################################################################################################
def __testOne(ID,oneCalc,engine='',calcType='',execPrefix=None,execPostfix=None): 
    """
    Stages the first the first step in a calculation workflow of a calc set to be submitted
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
    	  engine (str): executable that you are calling to run the calculations          
 	  execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)
	  execPostfix (str): commands to go after the executable when run
                            (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          calcType (str): used to pull the engine post processing executable to the calc dir
                          and to write the postprocessing input file if needed      

    Returns:
          None
          
    """

    if execPrefix == None:
        if AFLOWpi.prep.__ConfigSectionMap("run","execprefix") != '':
            execPrefix=AFLOWpi.prep.__ConfigSectionMap("run","execprefix")


        else:
            execPrefix=''

    if execPostfix == None:
        if AFLOWpi.prep.__ConfigSectionMap("run","execpostfix") != '':
            execPostfix = AFLOWpi.prep.__ConfigSectionMap("run","execpostfix")
        else:
            execPostfix=''


    logging.debug('entering __testOne') 
    try:
        global configFile
        configFile = oneCalc['_AFLOWPI_CONFIG_']
        AFLOWpi.prep.__forceGlobalConfigFile(configFile)

        if execPrefix == None:
            if AFLOWpi.prep.__ConfigSectionMap("run","execprefix") != '':
                execPrefix=AFLOWpi.prep.__ConfigSectionMap("run","execprefix")
            else:
                execPrefix=''

        if execPostfix == None:
            if AFLOWpi.prep.__ConfigSectionMap("run","execpostfix") != '':
                execPostfix = AFLOWpi.prep.__ConfigSectionMap("run","execpostfix")
            else:
                execPostfix=''
    except Exception,e:
        _fancy_error_log(e)

    try:
        engineDir  = AFLOWpi.prep.__ConfigSectionMap("prep",'enginedir')
        if os.path.isabs(engineDir) == False:
            configFileLocation = AFLOWpi.prep.__getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            engineDir =  os.path.join(configFileLocation, engineDir)
    except Exception,e:
        pass


    ri = ID+'.in'
    ro = ID+'.out' 

    rd = oneCalc['_AFLOWPI_FOLDER_']
    configFile = oneCalc['_AFLOWPI_CONFIG_']



    logging.debug('exiting __testOne')


#############################################################################################################

def __submitJob(ID,oneCalc,__submitNodeName__,sajOverride=False,forceOneJob=False):
    """
    Submits a step of a calculation's pipeline to be run
    
    Arguments:          
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          sajOverride (bool): Overrides stepsasjobs=False in the config file used
          forceOneJob (bool): Overrides stepsasjobs=True in the config file used
    Returns:
          None
          
    """

    logging.debug('entering __submitJob')
    folder = oneCalc['_AFLOWPI_FOLDER_']
    print 
#    try:
#        ID = __get_index_from_pp_step(oneCalc,ID)
#    except:
#        pass
    try:
        clusterType = AFLOWpi.prep.__ConfigSectionMap("cluster","type").upper()
        stepsAsJobs = AFLOWpi.prep.__ConfigSectionMap("cluster","stepsasjobs").lower()
        clusterType=clusterType.upper()

        if len(clusterType)==0:
            clusterType='None'
        elif clusterType==None:
            clusterType='None'
    except Exception,e:
        _fancy_error_log(e)


    """
    sets the cluster type to none if:
    1. the clusterType is PBS,
    2. the user has the flag stepsasjobs=false in their config file
    3. the job is being submitted from somewhere other than the submission node (i.e. one of the compute nodes)

    This will cause the first job (step in the calc chain) to be submitted to the cluster and all subsequent
    jobs will be run on the compute nodes as if it was the local machine. (so the calc chain is just one job
    on the cluster)
    """

    if socket.gethostname() == __submitNodeName__ and stepsAsJobs=='false':
        print 'Submitting steps as one cluster job per calculation chain.'
        logging.info('Submitting steps as one cluster job per calculation chain.')
        
    if (socket.gethostname() != __submitNodeName__ and stepsAsJobs.lower()!='true') or forceOneJob==True:
        if sajOverride == False:
            clusterType=''

###########################################################################################

    if ID == oneCalc['_AFLOWPI_PREFIX_'][1:]:
        pass
    else:
        ID='_'.join(ID.split('_')[:3])

#    global execFileString
    AFLOWpi.prep.__from_local_scratch(oneCalc,ID)

    try:
        #for calc subset...
        if oneCalc['_AFLOWPI_FOLDER_']==__main__.oneCalc['_AFLOWPI_FOLDER_']:
            #submit the cluster submission 
            submit_ID=__main__.ID
        else:
            submit_ID=ID
    except Exception,e:
#        AFLOWpi.run._fancy_error_log(e)
        submit_ID=ID

    if clusterType.upper() in ['PBS','SLURM']:
        '''transfer the .save folder from local scratch if the option is being used'''
#        AFLOWpi.prep.__from_local_scratch(oneCalc,ID,ext_list=['save','occup'])

        try:
            try:
                '''set parent processID to 0 so they'll be different when the script starts'''

                submitCommand='qsub'
                if clusterType.upper()=='SLURM':
                    submitCommand='sbatch'
                
                

            except Exception,e:
                _fancy_error_log(e)
#            qsubFilename = oneCalc['__qsubFileName__']
            
            
            qsubFilename = os.path.abspath(os.path.join(folder,'_'+submit_ID+'.qsub'))

            calcName=AFLOWpi.run.__get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])

            queue = ''
            if AFLOWpi.prep.__ConfigSectionMap("cluster","queue") !='':
                queue = '-q %s ' %  AFLOWpi.prep.__ConfigSectionMap("cluster","queue")
        except Exception,e:
            _fancy_error_log(e)
        try:
            
            command  = "ssh -o StrictHostKeyChecking=no %s '%s %s -N %s %s' " % (__submitNodeName__,submitCommand,queue,calcName,qsubFilename)
            job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True)
            job.communicate()               

            logging.info('submitted %s to the queue' % submit_ID)
        except Exception,e:
            try:
                command = "ssh -o StrictHostKeyChecking=no %s '%s %s' " % (__submitNodeName__,submitCommand,qsubFilename)
                job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True)
                job.communicate()               

            except Exception,e:
                _fancy_error_log(e)
                try:
                    baseID = ID.split('_')[0]
                    baseID = ID
                    qsubFilename = os.path.abspath(os.path.join(folder,'_'+baseID+'.qsub'))
                    os.system("ssh -o StrictHostKeyChecking=no %s '%s %s -N %s %s' " % (__submitNodeName__,submitCommand,queue,calcName,qsubFilename))
                except Exception,e:
                    print e
                    _fancy_error_log(e)

    else:
        try:
#            globals()['execFileString'] = os.path.abspath(os.path.join(folder,'_'+ID+'.py'))
            execFileString  = os.path.abspath(os.path.join(folder,'_'+ID+'.py'))
            

            execfile(execFileString)

        except KeyboardInterrupt:
            print 'Got Keyboard Exit signal, exiting'
            logging.debug('Got Keyboard Exit signal, exiting')
            sys.exit(0)
        except Exception,e:
            _fancy_error_log(e)

    logging.debug('exiting __submitJob')



##############################################################################################
import __main__
def submitFirstCalcs__(calcs):
    """
    Submits the first step of a calculation's pipeline
    
    Arguments:
          calcs (dict): Dictionary of dictionaries of calculations          

    Keyword Arguments:
          None

    Returns:
          None
          
    """


    logging.debug('entering __submitFirstCalcs')
    if '__submitNodeName__' not in globals().keys():
        global __submitNodeName__
        __submitNodeName__ = socket.gethostname()
    try:
        __submitNodeName__=__main__.__submitNodeName__
    except:
        pass
    logging.debug('sending from %s' % __submitNodeName__)
    
#    if '__submit__flag__' not in globals().keys():
#        logging.debug('exiting __submitFirstCalcs without submitting')
#        return
    if AFLOWpi.prep.__ConfigSectionMap('cluster','type') != '':
        print 'Submitting Jobs'
    for ID,oneCalc in calcs.iteritems():
        __submitJob(ID,oneCalc,__submitNodeName__)
    logging.debug('exiting __submitFirstCalcs')        


def submit():
    """
    sets global __submit__flag__ so calculations will start when user script completes
    
    Arguments:
          None

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    globals()['__submit__flag__']=True

#####################################################################################################################
def __restartPW(oneCalc,ID,a,__submitNodeName__,):
    """
    If pw.x ends becuause it hit max_seconds. set the input to restart the calculation
    and resubmit the job
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          a (float): Time that has passed since calculation has started
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    #mod input so it will not start from scratch next submission            
    logging.debug('Entering AFLOWpi.run.__restartPW')
    #    limiter=0.9
    try:
        walltime = AFLOWpi.run.__grabWalltime(oneCalc,ID)
    except:
        return 
#    try:
#        inFileString = AFLOWpi.retr.__writeInputFromOutputString(oneCalc,ID)
#    except:
    with open('%s.in' % ID,'r') as inputFileObj:
        inFileString= inputFileObj.read()
    inputDict = AFLOWpi.retr.__splitInput(inFileString)
    #bring wfc,hub, or any other files needed from local scratch if need be
    walltimeSec=90000000000
    try:
        walltimeSec = int(inputDict["&control"]["max_seconds"])
    except:
        '''if we can't find max_seconds in the input we just exit.'''
        logging.debug('Not pwscf. Exiting AFLOWpi.run.__restartPW')
        return
#        sys.exit(0)
    try:
        """see if calc reached it's end."""
        if (time.time()-a)>walltimeSec:
            inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
            """if it did then set to restart the calc"""

            inputDict["&control"]["restart_mode"] = "'restart'"
            ''' if the calc never had the time to start then just restart it from scratch'''

            inputDict["&control"]["max_seconds"] = int(walltimeSec)

            restartInput = AFLOWpi.retr.__joinInput(inputDict)
            '''write the new input file with the restarting'''
            with open('%s.in' % ID,"w") as inputfile:
                inputfile.write(restartInput)                  
            '''save to _<ID>.oneCalc'''
            oneCalc['_AFLOWPI_INPUT_']=restartInput
            AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
            """now resubmit"""
            resubmitID = oneCalc['__qsubFileName__'].split('/')[-1].split('.')[0][1:]
            #pull from local scratch before restarting

            #resubmit
            '''if we're at within 15 seconds of the end of walltime. don't submit and let restartScript do it'''

            oneCalc['__status__']['Restart']+=1

#            AFLOWpi.prep.__from_local_scratch(oneCalc,ID)
            AFLOWpi.run.__submitJob(resubmitID,oneCalc,__submitNodeName__,sajOverride=True)                             
            logging.debug('Restart Job Submitted. Exiting AFLOWpi.run.__restartPW')
            time.sleep(10)
            '''and exit cleanly.'''
            sys.exit(0)

        else:
            logging.debug('Job completed in time. No need for restart. Exiting AFLOWpi.run.__restartPW')

    except Exception,e:
        _fancy_error_log(e)


def __restartGIPAW(oneCalc,ID,a,__submitNodeName__,):
    """
    If pw.x ends becuause it hit max_seconds. set the input to restart the calculation
    and resubmit the job
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          a (float): Time that has passed since calculation has started
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    with open('%s.in' % ID,'r') as inputFileObj:
        inFileString= inputFileObj.read()
    inputDict = AFLOWpi.retr.__splitInput(inFileString)
    try:
        walltimeSec = int(inputDict["&control"]["max_seconds"])
    except:
        '''if we can't find max_seconds in the input we just exit.'''
        sys.exit(0)
    try:
        """see if calc reached it's end."""
        if (time.time()-a)>walltimeSec:
            inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
            """if it did then set to restart the calc"""

            inputDict["&inputgipaw"]["restart_mode"] = "'restart'"
            ''' if the calc never had the time to start then 
                just restart it from scratch'''
            try:
                if inputDict["&inputgipaw"]["max_seconds"]=='10':
                    inputDict["&inpitgipaw"]["restart_mode"] = "'from_scratch'"
            except:
                pass
            inputDict["&control"]["max_seconds"] = int(walltimeSec)
            restartInput = AFLOWpi.retr.__joinInput(inputDict)
            '''write the new input file with the restarting'''
            with open('%s.in' % ID,"w") as inputfile:
                inputfile.write(restartInput)                  
            '''save to _<ID>.oneCalc'''
            oneCalc['_AFLOWPI_INPUT_']=restartInput
            AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
            """now resubmit"""
            resubmitID = oneCalc['__qsubFileName__'].split('/')[-1].split('.')[0][1:]

            oneCalc['__status__']['Restart']+=1
            AFLOWpi.prep.__saveOneCalc(oneCalc,ID)

            __submitJob(resubmitID,oneCalc,__submitNodeName__,sajOverride=True)                 
            """give it a bit of time to get the resubmission through"""
#            time.sleep(10)
            '''and exit cleanly.'''
            sys.exit(0)
    except Exception,e:
        _fancy_error_log(e)



################################################################################################################
def __getWalltime(oneCalc,ID):
    """
    Get the walltime requested from the cluster submission script
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          walltime (int) amount of time requested
          
    """    

    cluster_type=AFLOWpi.prep.__ConfigSectionMap("cluster","type").upper()

    try:
        with open(oneCalc['__qsubFileName__'],"r") as qsubfileObj:
            qsubFileString = qsubfileObj.read()
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        return 50000000

    if cluster_type=='PBS':
        walltime = re.findall("walltime\s*=\s*([0-9:]*)",qsubFileString)[-1]

    if cluster_type=='SLURM':
        walltime = re.findall("--time\s*=\s*([0-9:]*)",qsubFileString)[-1]

    splitW = [int(x) for x in walltime.strip().split(":")]
    if len(splitW)==4:
        walltime=splitW[0]*24*3600+splitW[1]*3600+splitW[2]*60+splitW[3]
    elif len(splitW)==3:
        walltime=splitW[0]*3600+splitW[1]*60+splitW[2]
    elif len(splitW)==2:
        walltime=splitW[0]*60+splitW[1]
    elif len(splitW)==1:
        walltime=splitW[0]

    return walltime


def __generic_restart_check(oneCalc,ID,__submitNodeName__):
    """
    Checks to see if walltime limit is almost up. If it is within the buffer of time
    The job is resubmitted.
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          None
          
    """


    try:
        walltime,startTime=__grabWalltime(oneCalc,ID)    
    except:
        return

    percent_timer=0.90
    from_config = AFLOWpi.prep.__ConfigSectionMap("cluster","restart_buffer").lower()
    if from_config!='':
        percent_timer = float(from_config)
#    if percent_timer<1.0:
#        percent_timer=1.0-percent_timer

    
    if percent_timer<1.0:
        num_secs=int(float(walltime)*percent_timer)
    else:
        num_secs=int(float(walltime)-percent_timer)

    try:
        timediff=time.time()-startTime
    except Exception,e:
        timediff=0

    runTime = num_secs
    logging.debug('TIMEDIFF = %s'%timediff)      
    logging.debug('RUNTIME = %s'%runTime)
    logging.debug('num_secs = %s'%num_secs)
    logging.debug('percent_timer = %s'%percent_timer)
    logging.debug('STARTTIME = %s'%startTime)
    logging.debug('WALLTIME = %s'%walltime)

    '''if local scratch is being used don't run any steps if we're over 90% walltime'''
    if num_secs<(time.time()-startTime):
            '''set runTime to 0.0 make sure it trips the next if statement'''
            runTime=0.0
            

    if (runTime-60.0)<0.0:
        ID = __get_index_from_pp_step(oneCalc,ID)
        '''if there's not enough time to do anything so resubmit before anything runs'''
        '''make sure the scratch is copied back'''
        logging.info('job only has 60 seconds left. attempting to copy from local scratch if needed and resubmit job')
        '''force a new job'''
        AFLOWpi.run.__submitJob(ID,oneCalc,__submitNodeName__,sajOverride=True)                 
        logging.debug('Not enough time to start next step. Restarting')
        time.sleep(5)
        sys.exit(0)        



def __setupRestartPW(oneCalc,ID,__submitNodeName__):
    """
    Sets up pw.x input with max_seconds to it will cleanly exit if it
    doesn't finish within the time limit of the cluster submission 
    walltime requested.
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          
    """




    try:
        walltime,startTime=__grabWalltime(oneCalc,ID)
    except:
        '''if no walltime log in oneCalc (running local mode)'''
        return oneCalc,ID

    '''set buffer to time before walltime to stop the job'''
    percent_timer=0.90
    from_config = AFLOWpi.prep.__ConfigSectionMap("cluster","restart_buffer").lower()
    if from_config!='':
        percent_timer = float(from_config)


        
    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"r") as inputfile:
        inputStr=inputfile.read()
    '''try to split input to include max time if it doesn't work stop and return'''
    try:
        inputDict = AFLOWpi.retr.__splitInput(inputStr)
        calc_type=inputDict["&control"]["calculation"].strip(""" '`" """)
    except:
        return oneCalc,ID
    '''ibef calc type was found (i.e. not post processing) try to include max_seconds'''
    typeList = ['bands','nscf','scf','vc-relax','relax']

    try:
        timediff=time.time()-startTime
    except Exception,e:
        timediff=0

    if calc_type in typeList:
        try:
            max_sec=(int(float(walltime)*percent_timer))
            if percent_timer>1.0:
                max_sec=(int(float(walltime)-percent_timer))
            if  timediff<(walltime*percent_timer):
                runTime = max_sec-timediff                
                '''if script has already been running is small take the difference between walltime 
                and time script has been running. leave 10% of walltime to do cleanup.'''

                if (runTime-60.0)<0.0 or timediff>max_sec:
                    '''if there's not enough time to do anything so resubmit before anything runs'''
                    '''make sure the scratch is copied back'''
                    logging.info('Not enough time to start non-post processing calc. copying files from local scratch if needed and restarting')

                    ID = __get_index_from_pp_step(oneCalc,ID)
                    '''force a new job'''
                    AFLOWpi.run.__submitJob(ID,oneCalc,__submitNodeName__,sajOverride=True)                 
                    logging.debug('Not enough time to start pwscf. Restarting')
                    time.sleep(10)
                    sys.exit(0)        
            
            inputDict["&control"]["max_seconds"] = "%s" % int(runTime)
            restartInput = AFLOWpi.retr.__joinInput(inputDict)
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"w") as inputfile:
                inputfile.write(restartInput)
            oneCalc['_AFLOWPI_INPUT_']=restartInput
            AFLOWpi.prep.__saveOneCalc(oneCalc,ID)

        except Exception,e:
            logging.debug('if option to restart doesnt exist in the input of %s this will be thrown' % ID)
            _fancy_error_log(e)

        return oneCalc,ID



def __setupRestartGIPAW(oneCalc,ID):
    """
    Sets up GIPAW input with max_seconds to it will cleanly exit if it
    doesn't finish within the time limit of the cluster submission 
    walltime requested.
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          
    """

    try:
        walltime,startTime=__grabWalltime(oneCalc,ID)
    except:
        return oneCalc,ID
    
    globals()['__GLOBAL_CLOCK__']=startTime



    percent_timer=0.90
    try:
        percent_timer = AFLOWpi.prep.__ConfigSectionMap("cluster","restarttimer").lower()
        percent_timer = float(percent_timer.strip())
    except:
        percent_timer=0.90




    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"r") as inputfile:
        inputStr=inputfile.read()
    '''try to split input to include max time if it doesn't work stop and return'''
    try:
        inputDict = AFLOWpi.retr.__splitInput(inputStr)
        calc_type=inputDict["&inputgipaw"]["calculation"].strip(""" '`" """)
    except:
        return oneCalc,ID
    '''ibef calc type was found (i.e. not post processing) try to include max_seconds'''
    typeList = ['emr','efg','g_vectors','hyperfine']
    if calc_type in typeList:
        try:
            timediff=time.time()-globals()['__GLOBAL_CLOCK__']

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            timediff=0
        try:
            runTime = int(float(walltime)*percent_timer-timediff)
            if  timediff<(walltime*percent_timer):
                '''if script has already been running is small take the difference between walltime 
                and time script has been running. leave 10% of walltime to do cleanup.'''

            else:
                '''if there's not enough time to do anything so resubmit before anything runs'''
                return oneCalc,ID,False
            inputDict["&inputgipaw"]["max_seconds"] = "%s" % runTime
            restartInput = AFLOWpi.retr.__joinInput(inputDict)
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"w") as inputfile:
                inputfile.write(restartInput)
            oneCalc['_AFLOWPI_INPUT_']=restartInput
        except Exception,e:
            logging.debug('if option to restart doesnt exist in the input of %s this will be thrown' % ID)
            _fancy_error_log(e)

        return oneCalc,ID




def __writeWalltimeLog(oneCalc,ID,walltimeDict):
    """
    Writes the walltimeDict to disk ( _<ID>.oneCalc )
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          walltimeDict (dict): saves the walltime log to oneCalc and then to disk

    Keyword Arguments:
          None

    Returns:
          None
          
    """
    try:
        oneCalc['__walltime_dict__']=walltimeDict
        AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
    except:
        pass

def __readWalltimeLog(oneCalc,ID):
    """
    reads walltime log from oneCalc
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          walltimeLog (dict): dictionary containing start time for AFLOWpi as well
                              as the requested walltime for the cluster job 
          
    """

    return oneCalc['__walltime_dict__']






import __main__
def __grabWalltime(oneCalc,ID):
    """
    find the start time of current step in chain and record that into the walltime file
    
    Arguments:          
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation 

    Keyword Arguments:
          None

    Returns:
          walltimeStart (int): walltime requested
          startTime (int): start time of the script
          
    """

    try:        
        output = __readWalltimeLog(oneCalc,ID)

#        AFLOWpi.run._fancy_error_log(e)

        AFLOWpi.run.__writeWalltimeLog(oneCalc,ID,output)


        '''grab the first walltime in case stepsasjobs==false for all jobs in chain'''

        walltimeStart=output['walltime']    
        startTime=output['start']

        logging.debug('startTime Test: %s' %startTime)
    except Exception,e:
        logging.debug(e)
        raise SystemExit
        return float(40000000),float(0)
    return float(walltimeStart),float(startTime)

def __setStartTime(oneCalc,ID):
    """
    If this python script is the first to run in for AFLOWpi in a cluster job then
    the start time will be recorded and stored in the oneCalc object with the key  
    "__walltime_dict__" and the values for the start time of the script and the 
    requested walltime. 
    
    Arguments:          
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation 

    Keyword Arguments:
          None

    Returns:
          walltimeStart (int): walltime requested
          startTime (int): start time of the script
          
    """
    
    if AFLOWpi.prep.__ConfigSectionMap("cluster","type").lower()=='':
        return oneCalc
    try:
        try:
            walltime_dict=oneCalc['__walltime_dict__']
        except:
            oneCalc['__walltime_dict__']={}
            walltime_dict={}
    

        if  '__GLOBAL_CLOCK__' not in globals().keys():
            '''only move from local scratch if this script is starting fresh'''
            AFLOWpi.prep.__to_local_scratch(oneCalc,ID)
            oneCalc['__walltime_dict__']={}
            
            globals()['__GLOBAL_CLOCK__']=__main__.__CLOCK__
            walltime_dict['walltime'] = AFLOWpi.run.__getWalltime(oneCalc,ID)         

            walltime_dict['start']=__main__.__CLOCK__
            AFLOWpi.run.__writeWalltimeLog(oneCalc,ID,walltime_dict)

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

        try: 
            oneCalc['__walltime_dict__']={}
            walltime_dict={}
            walltime_dict['walltime'] = 40000000000.0         
            walltime_dict['start']=0.0
            AFLOWpi.run.__writeWalltimeLog(oneCalc,ID,walltime_dict)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
        
    return oneCalc


def __restartScript(oneCalc,ID,PID):
    '''
    special script started when the _<ID>.py script starts. reads the walltime in the qsub file and 
    resubmits and then exits if the time ran gets within 1 minute of the walltime requested. 
    Used for when the script is on the post processing stages and may happen to get killed by
    the cluster daemon.
    '''

    '''keep trying to find oneCalc and ID from '''


    logging.debug('main processess running on PID: %s' % PID)
    percent_timer=0.90
    try:
        percent_timer = AFLOWpi.prep.__ConfigSectionMap("cluster","restarttimer").lower()
        percent_timer = float(percent_timer.strip())
    except:
        percent_timer=0.90

    walltime,startTime = __grabWalltime(oneCalc,ID)
    globals()['__GLOBLAL_CLOCK__']=startTime
    
    walltime  =  int(walltime)
#    __main__.__GLOBAL_CLOCK__=float(startTime)
#    globals()['__GLOBAL_CLOCK__']=startTime
    
    __submitNodeName__=__main__.__submitNodeName__

    stepsasjobs = AFLOWpi.prep.__ConfigSectionMap("cluster","stepsasjobs").lower()
    if stepsasjobs=='false' or stepsasjobs=='f':
        stepsasjobsBool=False

        logging.debug('stepsasjobs set to false in config file. reading in start time from oneCalc')

    currentTime=time.time()-startTime
    logging.debug('walltime in seconds: %s' % walltime)

    logging.debug('Starting monitoring thread for %s in %s' %(ID,oneCalc['_AFLOWPI_FOLDER_']))
    remainingTime = float(walltime)-float(currentTime)
    logging.debug('time since start for _<ID>.py = %s' %  currentTime)
    logging.debug('current time till ending _<ID>.py = %s' %  remainingTime)
    '''restart script and exit immediately if it's about to be killed by the cluster'''
    while walltime-currentTime > 5:
        '''update time to reflect current time'''        
        currentTime=time.time()-globals()['__GLOBAL_CLOCK__']
        '''
        if we've moved into another step with stepsasjobs==false then we don't 
        need the restartScript from this step anymore so we should exit the thread
        '''
        try:
            if globals()['__TEMP__INDEX__COUNTER__'] != oneCalc['__chain_index__']:
                return
        except Exception,e:
            logging.debug(e)
        time.sleep(1.0)

        try:
            '''check to see if main script is running..if not it's 
           not then something went wrong in the main script and this should also exit'''

            os.kill(PID, 0) 
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            logging.info('main script did not exit cleanly. killing __restartScript subprocess')
            print 'main script did not exit cleanly. killing __restartScript subprocess'
            sys.exit(0)


    logging.info('Job not completed in time. Restarting from that step.')

    oneCalc =AFLOWpi.prep.__loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],ID)
    try:
        oneCalc['__status__']['Restart']+=1
        AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
        AFLOWpi.run.__submitJob(ID,oneCalc,__submitNodeName__,sajOverride=True)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
    logging.info('Resubmission successful. Waiting 5 seconds then ending main script')
    '''after the job is submitted wait ten seconds to exit the script'''
    time.sleep(60)
    '''and exit gracefully'''
    logging.info('telling main process to stop')
    os.kill(PID,signal.SIGUSR1) 
    sys.exit(0)


##################################################################################################################
def __convert_fortran_double(fort_double_string,string_output=False):
    """
    find the start time of current step in chain and record that into the walltime file
    
    Arguments:
          fort_double_string (str): Fortran double in string form

    Keyword Arguments:
          string_output (bool): Output a string or a float

    Returns:
          fort_double_string (float): the fotran double in float form or string
          
    """

    fort_double_string=fort_double_string.lower()
    fort_double_list=fort_double_string.split('d')
    if len(fort_double_list)==1:
        fort_double_list=fort_double_string.split('f')
        if len(fort_double_list)==1:
            print 'could not convert fortran float string to numpy float'
            logging.info('could not convert fortran float string to numpy float')
            return fort_double_string
    if len(fort_double_list)>2:
        print 'could not convert fortran float string to numpy float'
        logging.info('could not convert fortran float string to numpy float')
        return fort_double_string
    try:

        mult=10**int(fort_double_list[1])
        base=float(fort_double_list[0])
        retr_float=base*mult
        if string_output==True:
            return repr(retr_float)
        else:
            return retr_float

    except Exception,e:
        print e
        print 'could not convert fortran float string to numpy float'
        logging.info('could not convert fortran float string to numpy float')
        return fort_double_string

# def __lower_mixing(oneCalc,ID):


#     inputDict=AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
#     try:
#         if '&electrons' not in inputDict.keys():
#             inputDict['&electrons']=collections.OrderedDict()
#         old_beta='0.7d0'
#         if 'mixing_beta' in inputDict['&electrons'].keys():
#             old_beta=inputDict['&electrons']['mixing_beta']

#         old_beta_float=__convert_fortran_double(old_beta)
#         new_beta_float=old_beta_float-0.2
#         if new_beta_float>0.09:
#             inputDict['&electrons']['mixing_beta']=str(new_beta_float)
#             inputString=AFLOWpi.retr.__splitInput(inputDict)
#             oneCalc['_AFLOWPI_INPUT_']=inputString

#             outFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in')
#             with open(outFileName,'w') as outFile:
#                 outFile.write(inputString)
#     except:
#         pass
#     AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
    
#     return oneCalc,ID
    

def __get_index_from_pp_step(oneCalc,ID):
    """
    DEFUNCT. CANDIDATE FOR REMOVAL
    
    Arguments: 
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation         

    Keyword Arguments:
          None  

    Returns:
          ID (str): Returns the input ID
          
    """

    return ID
    try:
        index = int(oneCalc['__chain_index__'])
        if index==1 or index==0:
            return oneCalc['_AFLOWPI_PREFIX_'][1:]
        else:
            return '%s_%02d'%(oneCalc['_AFLOWPI_PREFIX_'][1:],index)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        return ID


def __PW_bands_fix(oneCalc,ID):
    """
    Accounts for the occational problem of the bands.x output being 
    improperly formatted for AFLOWpi to parse later.
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          None

    Returns:
          None 
          
    """


    rd = os.path.abspath(oneCalc['_AFLOWPI_FOLDER_'])

    try:
        """fix bands.out where you have large numbers butting up against each other ex.) -113.345-112.567 
        and split the two numbers up for processing with band_plot.x or else you'll get an error"""
        with open(os.path.join(rd,'%s_band_data.out'%ID),'r') as fileinput:
            input = fileinput.read()
            with open(os.path.join(rd,'%s_band_data.out'%ID),'w') as fileoutput:
                correctedText = re.sub(r'(\d)([-][\d])',r'\1 \2',input)
                fileoutput.write(correctedText)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)


    try:
        ##### make the bands.xmgr now that band_plot.x is gone
        with open(os.path.join(rd,'%s_band_data.out'%ID),'r') as bands_out_file:
            bands_out_lines = bands_out_file.readlines()

        indent_re=re.compile(r'      ')
        total=0.0
        path_vals=[]
        k_val=[]
        total_list=[]
        for i in bands_out_lines:
            if indent_re.match(i)!=None:
                try:
                    dist=numpy.asarray([float(x) for x in i.split()])-total_path_length
                    path_vals.append(k_val)
                    total_list.append(total)            
                    total+=numpy.linalg.norm(dist)

                    k_val=[]
                except Exception,e:
                    pass
                total_path_length=numpy.asarray([float(x) for x in i.split()])
            else:
                try:
                    k_val.extend(numpy.asarray([float(x) for x in i.split()]))
                except Exception,e:
                    pass

        path_vals.append(k_val)
        total_list.append(total)
        path_vals=numpy.array(path_vals)
        output_string=''
        for item in path_vals.T:
            print_list=[]
            for j in range(len(item)):
                output_string+='%12.7f %12.5f\n'%(total_list[j],item[j])
            output_string+='\n'

        with open(os.path.join(rd,'%s_bands.xmgr'%ID),'w') as bands_out_file:
            bands_out_file.write(output_string)

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)




def __vcrelax_error_restart(ID,oneCalc,__submitNodeName__):
    """
    Restarts if vc-relax calculation fails

    BROKEN/NOT NEEDED POSSIBLY DELETE
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          None
          
    """

    error_regex_list=[]


    submit_ID = __get_index_from_pp_step(oneCalc,ID)
    try:        
        if AFLOWpi.prep.__check_restart(oneCalc,ID):
            oneCalc,submit_ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','restart_mode',"'from_scratch'")
    except:
        pass
    with open('%s.out'% ID,'r') as outFileObj:
        outFileString = outFileObj.read()
        for re_error in  error_regex_list:
            if len(re_error.findall(outFileString))!=0:
                try:
                    pass
                except:
                    pass
                try:
                        oneCalcBase =AFLOWpi.prep.__loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],submit_ID)
                        oneCalcBase['__runList__']=[]
                        oneCalc['__runList__']=[]
                        AFLOWpi.prep.__saveOneCalc(oneCalcBase,submit_ID)
                except:
                    pass

                AFLOWpi.run.__submitJob(submit_ID,oneCalc,__submitNodeName__,forceOneJob=True)
                

def __io_error_restart(ID,oneCalc,__submitNodeName__):
    """
    Restarts if I/O error encountered

    BROKEN/NOT NEEDED POSSIBLY DELETE
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          None
          
    """


    try:
        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','restart_mode',"'from_scratch'")
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    error_regex_list=[]
    error_regex_list.append(re.compile(r'Error in routine diropn'))
    error_regex_list.append(re.compile(r'Error in routine davcio'))
    error_regex_list.append(re.compile(r'Error in routine pw_readfile'))
    error_regex_list.append(re.compile(r'Error in routine seqopn'))
    error_regex_list.append(re.compile(r'kpt grid not Monkhorst-Pack'))
    submit_ID = __get_index_from_pp_step(oneCalc,ID)
    with open('%s.out'% ID,'r') as outFileObj:
        outFileString = outFileObj.read()
        for re_error in  error_regex_list:
            if len(re_error.findall(outFileString))!=0:
                if re_error.findall(outFileString)=='Error in routine pw_readfile' or re_error.findall(outFileString)=='kpt grid not Monkhorst-Pack':
                    try:
                        oneCalcBase =AFLOWpi.prep.__loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],submit_ID)
                        oneCalcBase['__runList__']=[]
                        oneCalc['__runList__']=[]
                        AFLOWpi.prep.__saveOneCalc(oneCalcBase,submit_ID)
                    except:
                        pass
                try:

                    logging.warning("%s: Trying again..."%re_error.findall(outFileString)[-1])
                except:
                    pass
                
                try:
                    oneCalcBase =AFLOWpi.prep.__loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],submit_ID)

                        
                    restart_num = oneCalcBase['__status__']['Restart']
                    oneCalcBase['__status__']['Restart']=restart_num+1
                    oneCalc['__status__']['Restart']=restart_num+1
                    if restart_num>10:
                        oneCalcBase["__execCounter__"]=0
                        oneCalc["__execCounter__"]=0
                    if restart_num>15:
                        sys.exit(0)

                    AFLOWpi.prep.__saveOneCalc(oneCalcBase,submit_ID)
                except:
                    pass

                AFLOWpi.run.__submitJob(submit_ID,oneCalc,__submitNodeName__,forceOneJob=True)
                sys.exit(0)

##################################################################################################################

def __qe__pre_run(oneCalc,ID,calcType,__submitNodeName__,engine):
    """
    Performs the pre-run checks and restarts the cluster job if needed. 
    On local mode this gets skipped
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation
          calcType (str): type of calculation
          __submitNodeName__ (str): String of hostname that cluster jobs should be submitted from

    Keyword Arguments:
          None

    Returns:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation          
          
    """

    clusterType = AFLOWpi.prep.__ConfigSectionMap("cluster","type")
    home = os.path.abspath(os.curdir)
    rd ='./'

    if calcType=='scf':    
        oneCalc,ID = AFLOWpi.prep.__setup_local_scratch(oneCalc,ID)

    try:
        if clusterType.upper() in ['PBS','SLURM']:


            #         #if it's scf/nscf setup local scratch
            #         #if it's a restart files will be transferred to nodes

            # else:
            #         try:

            #             #if it's post processing pull the files from the local scratch
            #             pass

            #         except Exception,e:
            #             AFLOWpi.run._fancy_error_log(e)

            '''if there's 60 seconds or less left exit immediately and don't start a new job'''

            AFLOWpi.run.__generic_restart_check(oneCalc,ID,__submitNodeName__)
            if calcType in ['scf']:    
                oneCalc,ID = AFLOWpi.run.__setupRestartPW(oneCalc,ID,__submitNodeName__)



            if calcType in ['emr','nmr','hyperfine','gvectors']:
                oneCalc,ID,restartBool = __setupRestartGIPAW(oneCalc,ID)




    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)


    try:
        if calcType!='scf' and calcType!=None and calcType!='custom':
                with open(os.path.join(rd,"%s_%s.in") % (ID,calcType),'w+') as PPInput:
                    PPInput.write(__makeInput(oneCalc,engine,calcType,ID=ID))

    except Exception,e:
        _fancy_error_log(e)



    return oneCalc,ID

###############################################################################################################

def __oneRun(__submitNodeName__,oneCalc,ID,execPrefix='',execPostfix='',engine='espresso',calcType=None,executable=None,execPath=None,nextCalc=None,nextConf=None,):
    """
    Run a single calculation in the dictionary with a specific engine
    
    Arguments:

	  oneCalc (dict): dictionary of the calculation
          ID (str): Identifying hash of the calculation
    Keyword Arguments:
	  engine (str): executable that you are calling to run the calculations          
 	  execPrefix (str): commands to go before the executable when run 
                            (ex. mpiexec nice -n 19 <executable>) (default = None)
	  execPostfix (str): commands to go after the executable when run
                            (ex. <execPrefix> <executable> -ndiag 12 -nimage 2) (default = None)
          calcType (str): used to pull the engine post processing executable to the calc dir
                          and to write the postprocessing input file if needed
          executable (str): if the executable has already been copied to the calc's directory
                            then use this to run the input <ID>.in
          execPath (str): path of the executable if needed
          nextCalc (str): DEFUNCT
          nextConf (str): DEFUNCT

          
    Returns:
          None
          
    """

    logging.debug('entering __oneRun')


    os.chdir(os.path.abspath(oneCalc['_AFLOWPI_FOLDER_']))

    AFLOWpi.prep.__forceGlobalConfigFile(oneCalc['_AFLOWPI_CONFIG_'])
    clusterType = AFLOWpi.prep.__ConfigSectionMap("cluster","type")

    try:
        '''save the status of the calcs'''
        if '__status__' not in oneCalc.keys():
            oneCalc['__status__']=collections.OrderedDict()
        oneCalc['__status__']['Start']=True
        oneCalc['__status__']['Complete']=False
        AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
    # '''set the completed status to false in the case of a looping set of calcs'''
    # 





    home = os.path.abspath(os.curdir)
    rd ='./'

    '''do the setup'''

    oneCalc,ID = AFLOWpi.run.__qe__pre_run(oneCalc,ID,calcType,__submitNodeName__,engine)


#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
###CONSIDER MOVING OUTSIDE
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
        ########################################################################################################





       ########################################################################################################


       ########################################################################################################



#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
###CONSIDER MOVING OUTSIDE
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
    try:
        if calcType!='scf' and calcType!=None and calcType!='custom':
            ri = ID+'_'+calcType+'.in'
            ro = ID+'_'+calcType+'.out'
        else:
            ri = ID+'.in'
            ro = ID+'.out'




        if execPath==None:
            execPath=''
            executable = __getExecutable(engine,calcType)

            executable = os.path.join('./',executable)
        else:
            execPath = execPath
            executable=''

        if os.path.exists(os.path.join(rd,engine)):
            command='%s %s%s %s < %s  >  %s' % (execPrefix,execPath,executable,execPostfix,ri,ro)
            logging.error(command)
        else:
            command='%s %s%s %s < %s  >  %s' % (execPrefix,execPath,executable,execPostfix,ri,ro)

    except Exception,e:
        _fancy_error_log(e)
###############################################################################################################      
    try:	
        logging.info("starting %s in %s" % (command, home))
        print "starting %s in %s" % (command, home)

        if clusterType.upper() in ['PBS','SLURM']:
            job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True,cwd=rd)
            a = time.time()

            job.communicate()               

                

#                typeList = ['band]
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
###CONSIDER MOVING OUTSIDE
#################################################################################################################      
#################################################################################################################      
#################################################################################################################     
            if clusterType.upper() in ['PBS','SLURM']:
                AFLOWpi.run.__restartPW(oneCalc,ID,a,__submitNodeName__) 

            if job.returncode==0:
                try:
#                    print job.returncode
                    if calcType in ['scf']:
                        pass
                    elif calcType in ['emr','nmr','hyperfine','gvectors']:
                        __restartGIPAW(oneCalc,ID,a,__submitNodeName__)
                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)
               # __io_error_restart(ID,oneCalc,__submitNodeName__)



            if job.returncode==1 or job.returncode==255:
                logging.warning('%s in %s returned exit code: %s'%(ID,oneCalc['_AFLOWPI_FOLDER_'],job.returncode))

                try:
                    if oneCalc['__status__']['Error']!='None':
                        oneCalc['__status__']['Error']+='Exited Status: %s' % 1
                    else:
                        oneCalc['__status__']['Error']='Exit Status: %s' % 1
                except:
                    oneCalc['__status__']['Error']='Exited Status %s' % 1 
                AFLOWpi.prep.__saveOneCalc(oneCalc,ID)                                


                sys.exit(0)

            elif job.returncode==2:
                oneCalc['__status__']['Error']='FAILED TO CONVERGE %s' % job.returncode
                AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
                sys.exit(0)

            elif job.returncode > 2 and job.returncode <= 127:				
                #empty text file is created in the AFLOWpi folder
                #is created and the keys of the failed to converge
                #calculations is written in the file
                oneCalc['__status__']['Error']='QE ERR %s' % job.returncode
                AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
                sys.exit(0)
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
###CONSIDER MOVING OUTSIDE
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      

        else:
            os.system(command)
       ###############################################################################################################
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
###CONSIDER MOVING OUTSIDE
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
        logging.info("finished %s in %s" % (command, rd))
        print "finished %s in %s" % (command, rd)

        if AFLOWpi.prep.__ConfigSectionMap('prep','savedir') != '':
            try:

                AFLOWpi.retr.__moveToSavedir(os.path.join(rd,ro))
                AFLOWpi.retr.__moveToSavedir(os.path.join(rd,ri))
                logs = glob.glob('%s/*%s*.log' % (rd,oneCalc['_AFLOWPI_PREFIX_'][1:]))
                for log in logs:
                    AFLOWpi.retr.__moveToSavedir(log)
                pdf = glob.glob('%s/*%s*.pdf' % (rd,oneCalc['_AFLOWPI_PREFIX_'][1:]))
            except Exception,e:
                _fancy_error_log(e)

        logging.debug('STARTING DATA')
        if oneCalc['_AFLOWPI_PREFIX_'][1:]==ID:
            logging.debug('RECORDING DATA')
            try:

                AFLOWpi.retr.getCellOutput(oneCalc,ID)

            except Exception,e:
                pass
       ###############################################################################################################
    except KeyboardInterrupt:
        print 'Got Keyboard Exit signal, exiting'
        logging.debug('Got Keyboard Exit signal, exiting')
        sys.exit(0)
    except OSError as e:
            print "ERROR! "+str(e)
            logging.error("ERROR! "+str(e))
    except ValueError as e:
            print "ERROR! "+str(e)
            logging.error("ERROR! "+str(e))
    except subprocess.CalledProcessError,e:
            if e.returncode==1:
                print "exit code==1"


        ##############################################################################################################

#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
###CONSIDER MOVING OUTSIDE
#################################################################################################################      
#################################################################################################################      
#################################################################################################################      
    os.chdir(home)		



    oneCalc['__status__']['Complete']=True
    AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
    logging.debug('exiting __oneRun')





def __onePrep(oneCalc,ID,execPrefix=None,execPostfix=None,engine='espresso',calcType='',alt_ID=None):
    """
    Prepares one calculation run an engine executable
    
    Arguments:
          oneCalc (dict): dictionary of one of the calculations
          ID (str): ID of calculation

    Keyword Arguments:
          engine (str): executable that you are calling to run the calculations (default='pw.x')          
          execPrefix (str): commands to go before the executable when run (ex. mpiexec nice -n 19 <executable>) 
                            (default = None)
          execPostfix (str): commands to go after the executable when run (ex. <execPrefix> <executable> -npool 12 )
                             (default = None)
          calcType (str): used to prep for a particular type of calc (i.e. 'scf','pdos'..etc) 

    Returns:
          None
          
    """

    logging.debug('entering __onePrep')

    if execPrefix == None or execPrefix ==  '':
            execPrefix=AFLOWpi.prep.__ConfigSectionMap("run","execprefix")
    else:
        execPrefix=''
    if execPostfix == None or execPostfix ==  '':

        execPostfix = AFLOWpi.prep.__ConfigSectionMap("run","execpostfix")
    else:
        execPostfix=''

    folder = oneCalc['_AFLOWPI_FOLDER_']
    prefix = oneCalc['_AFLOWPI_PREFIX_']


    """run projwfc.x on all the folders"""
    home = os.path.abspath(os.curdir)
    rd = os.path.abspath(folder)
    """dos.out is the output from dos.x as designated within dos.in"""


    try:
        if alt_ID==None:
            write_ID='ID'
        else:
            write_ID="'%s'"%alt_ID

        execString=''
        execString+='''if oneCalc['__execCounter__']<=%s:
        AFLOWpi.run.__oneRun''' % oneCalc['__execCounterBkgrd__']

        execString+='''(__submitNodeName__,oneCalc,%s,execPrefix='%s',execPostfix='%s',engine='%s',calcType='%s')
        oneCalc['__execCounter__']+=1
        AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
        '''  % (write_ID,execPrefix,execPostfix,engine,calcType)

        oneCalc['__execCounterBkgrd__']+=1

        execFileString = os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py'))
        logging.debug('ID: %s' % ID)


        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN',execString)
        AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN',"'%s'" % ID)


    except Exception,e:
        print e
        AFLOWpi.run._fancy_error_log(e)


    logging.debug('exiting __onePrep')
################################################################################################################

################################################################################################################
def __makeInput(oneCalc,engine,calcType,ID=''):
    """
    Writes the input file for postprocessing calculations
    
    Arguments:
          oneCalc (dict): Dictionary of one of the calculations          
          engine (str): Particular engine executable being used for the postprocessing step 
                        that you are calling to run the calculations (default='pw.x')          
          calcType (str): Type of PP calculation to be done
          
    Keyword Arguments:
          ID (str): ID hash for the calculation. Needed for the filename (<ID>_<calcType>.in)

    Returns:
          stringDict (str): Input string of the PP step
          
    """

    try:
        prefix = AFLOWpi.retr.__prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])
    except:
        prefix = oneCalc['_AFLOWPI_PREFIX_']

    temp_dir= AFLOWpi.prep.__get_tempdir()
    prefix=prefix.replace('"','') 
    prefix=prefix.replace("'","") 
    calcType=calcType.lower()

    if engine=='espresso':
        bandsInput=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID)
        nbnd=''
        if calcType!='scf':
            bandsInput = oneCalc['_AFLOWPI_INPUT_']
        else:
            nbnd='None'
        try:
            nbnd = re.findall(r'nbnd\s*=\s*([0-9]+)',bandsInput)[-1]
        except Exception,e:
            nbnd='30'

        stringDict={'dos':
                        """ &dos
       prefix='%s'
       outdir='%s'
       fildos='%s_dos.dat'

      /
""" % (prefix,temp_dir,ID),
                            'pdos':"""  &projwfc
       prefix='%s'
    !  DeltaE=0.03
!       Emax=20
!       Emin=-20
       outdir='%s'

    !    kresolveddos=.true.
       filpdos='%s'

      / 
 """ % (prefix,temp_dir,ID),

                    'bands':""" &bands
     prefix='%s'
     outdir='%s'
     filband='./%s_band_data.out'
 /
 """ %  (prefix,temp_dir,ID),


                    'scf':
                        'None',
                    'nmr':
                        '''&inputgipaw
  job = 'nmr'
  prefix='%s'

/
''' %prefix,
                    'g_tensor':
                        '''&inputgipaw
  job = 'g_tensor'
  prefix='%s'

/''' %prefix,

                    'emr':
                        '''&inputgipaw
  job = 'efg'
  prefix='%s'
                    
/

''' %prefix,
                    'hyperfine':
                        '''&inputgipaw
  job = 'hyperfine'
  prefix='%s'
  hfi_output_unit = 'MHz'
/
''' %prefix,


                            }

        

    return stringDict[calcType]
      
################################################################################################################        
################################################################################################################
 
        
            


                 

        

