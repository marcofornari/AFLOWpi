import AFLOWpi
import ast
import functools
import os
import subprocess
import glob
import re
import numpy
import collections

def _one_phon(oneCalc,ID):
    return glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/FD_PHONON'+'/displaced*.in')


def _phonon_band_path(oneCalc,ID,nk=400):
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

    dk=0.0001
    nk=10000
    path = AFLOWpi.retr._getPath(dk,oneCalc,ID=ID)

    '''scale the k points in the path list so they're as close to nk as possible'''
#    splitPath =  [[x.split()[3],int(x.split()[-1])] for x in  path.split('\n')[1:] if len(x)!=0]
    total =  [int(x.split()[3]) for x in  path.split('\n')[1:] if len(x)!=0]

    for entries in range(len(total)):
	    if total[entries]==0:
		    total[entries]+=1

    total=sum(total)
    scaleFactor = float(nk)/float(total)

    dk_prime = dk/scaleFactor

    path = AFLOWpi.retr._getPath(dk_prime,oneCalc,ID=ID)

    return path

def _pull_forces(oneCalc,ID):
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
        return
    force_split = [x.split()[6:] for x in force_block.split('\n') if len(x.strip())!=0]

    force_out_string=''
    for i in force_split:
        force_out_string+='%22.18s%22.18s%22.18s\n'%(i[0],i[1],i[2],)

    force_postfix='.'.join(ID.split('.')[1:])
    force_postfix=force_postfix.split('_')[0]

    with open(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'force.%s'%force_postfix),'w+') as out_file_obj:
        out_string = out_file_obj.write(force_out_string)




def _pp_phonon(__submitNodeName__,oneCalc,ID,LOTO=True,de=0.01,raman=True,field_strength=0.001,project_phDOS=True):
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

        

    if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
	    execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
    else:
        execPrefix=''
            
    #run fd_ifc
    AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_fd_ifc'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./fd_ifc.x' ) 

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

    epol_list=glob.glob("./FD_PHONON/epol.*")
    if len(epol_list)!=0:
      try:

        for field in field_list:
            try:

                AFLOWpi.run._gen_fd_input(oneCalc,ID,for_type=field,de=field_strength)
                AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_%s'%(ID,field),execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./fd_ef.x' )


            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)            

#######################################################################
        if raman==True or LOTO==True:
            epol_string=AFLOWpi.run._pull_eps_out(oneCalc,ID)
            born_charge_string=AFLOWpi.run._pull_born_out(oneCalc,ID)

            LOTO_REPLACE='''T

    %s
    %s

    ''' % (epol_string,born_charge_string)


            header_string=re.sub("F\n",LOTO_REPLACE,header_string)
########################################################################
      except Exception,e:
          AFLOWpi.run._fancy_error_log(e)

    else:
        pass
    
    ifc_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_ifc.fc'%ID)
    with open(ifc_file,'r') as ifc_fileobj:
        ifc_string=ifc_fileobj.read()

    for_matdyn  = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.fc' % ID)
    with open(for_matdyn,'w') as mat_dyn_in_fileobj:
        mat_dyn_in_fileobj.write(header_string)
        mat_dyn_in_fileobj.write(ifc_string)


    #run matdyn for bands
    AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_matdyn_phBand'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )

    #run matdyn for dos
    AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_matdyn_phDOS'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )


    if project_phDOS:
        try:
            AFLOWpi.run._project_phDOS(oneCalc,ID) 
            appdos_fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.eig.ap"%ID)
            AFLOWpi.run._sum_aproj_phDOS(appdos_fn,de=2.0)
        except Exception,e:
    	    AFLOWpi.run._fancy_error_log(e)
    #remove the fd.x output files in case we do another phonon calc
    try:
	    globpath=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'FD_PHONON')
	    phil = glob.glob(globpath+'/displaced*.in')
	    for i in phil:
		    os.remove(i)
    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)


def write_fdx_template(oneCalc,ID,nrx1=2,nrx2=2,nrx3=2,innx=2,de=0.01,atom_sym=True,disp_sym=True,proj_phDOS=True):
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

    temp_dir= AFLOWpi.prep._get_tempdir()

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
      file_out='./%s_ifc'
      noatsym        = %s
      nodispsym      = %s
      verbose=.true.
      hex=.false.
    /
     '''%(oneCalc['_AFLOWPI_PREFIX_'],nrx1,nrx2,nrx3,innx,de,ID,noatsym,nodispsym)

    FD_ifc_in_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_fd_ifc.in'%ID)
    with open(FD_ifc_in_path,'w') as fd_ifc_in_file:
        fd_ifc_in_file.write(fd_ifc_template)




    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

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


    ph_band_path=AFLOWpi.run._phonon_band_path(oneCalc,ID)



    matdyn_template='''
 &input
   asr='crystal',
%s
   flvec='%s.phBAND.modes'
   flfrc='%s.fc',
   flfrq='%s.phBAND', 
   fldyn='%s.phBAND.dyn',
   fleig='%s.phBAND.eig',
!   l1=%d
!   l2=%d
!   l3=%d
!   nosym=%s
   q_in_cryst_coord = .true. 
   fd=.true.
   na_ifc=.true.
   q_in_band_form=.true.,
   eigen_similarity=.true.

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


    #don't use sym q point reduction if projecting so that all weights are the same
    if proj_phDOS:
        for_proj_phDOS='.TRUE.'
    else:
        for_proj_phDOS='.FALSE.'

    matdyn_dos_template='''
 &input
   asr='crystal',
%s
   fd=.true.
   na_ifc=.true.
   fldos='%s.phdos'
   flfrc='%s.fc',
   flfrq='%s.phDOS', 
   fleig='%s.phDOS.eig'
   fldyn='%s.phDOS.dyn'
   flvec='%s.phDOS.modes'
   nosym=%s
!   l1=%d
!   l2=%d
!   l3=%d
   dos=.true.,
   eigen_similarity=.false.
%s   
 /
'''%(amass_str,ID,ID,ID,ID,ID,ID,for_proj_phDOS,nrx1,nrx2,nrx3,kgrid)
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


def prep_fd(__submitNodeName__,oneCalc,ID,nrx1=2,nrx2=2,nrx3=2,innx=2,de=0.01,atom_sym=True,disp_sym=True,proj_phDOS=True):

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
#    if ID in oneCalc["prev"]:
#        return

    #copy fd.x to the directories
    engineDir  = AFLOWpi.prep._ConfigSectionMap("prep",'engine_dir')
    # make sure the .save directory is copied back from local scratch
    # if that option is being used
    AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.save'])
    AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.occup'])

    for fd_exec in ['fd.x','fd_ifc.x','fd_ef.x','matdyn.x']:
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
            AFLOWpi.prep.totree(os.path.join(engineDir,'%s'%fd_exec),{ID:oneCalc},symlink=False)
        else:
            AFLOWpi.prep.totree(os.path.join(engineDir,'%s'%fd_exec),{ID:oneCalc},symlink=True)

    write_fdx_template(oneCalc,ID,nrx1=nrx1,nrx2=nrx2,nrx3=nrx3,innx=innx,de=de,atom_sym=atom_sym,disp_sym=disp_sym,proj_phDOS=proj_phDOS)
    


    AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_fd'%ID,execPrefix='',execPostfix='',engine='espresso',calcType='custom',execPath='./fd.x' )    


    ocd = oneCalc['_AFLOWPI_FOLDER_']
    globpath=os.path.join(ocd,'FD_PHONON')

    phil = glob.glob(globpath+'/displaced*.in')



    infile=    reduce_kpoints(oneCalc["_AFLOWPI_INPUT_"],[nrx1,nrx2,nrx3,])
    splitInput = AFLOWpi.retr._splitInput(infile)
    kpoints = splitInput['K_POINTS']['__content__']
    mod = splitInput['K_POINTS']['__modifier__']

    KPS="K_POINTS %s\n%s\n"%(mod,kpoints)


    for i in phil:
        fdfo = open(i,'r')
        fdfos = fdfo.read()
        fdfo.close()
        NC=clean_cell_params(fdfos)

        NC+="\n"+KPS
        fdfo = open(i,'w')
        fdfo.write(NC)
        fdfo.close()

def reduce_kpoints(inputfile,factor):
    try:    
                    splitInput = AFLOWpi.retr._splitInput(inputfile)
                    mod = splitInput['K_POINTS']['__modifier__'].upper()
                    mod=mod.strip('{}()')
                    if mod=='GAMMA':
                        inputfile = AFLOWpi.retr._joinInput(splitInput)
                    else:
                        splitInput=AFLOWpi.retr._splitInput(inputfile)
                        scfKPointString = splitInput['K_POINTS']['__content__']
			scfKPointSplit = [float(x) for x in scfKPointString.split()]


                            
                            

			for kpoint in range(len(scfKPointSplit)-3):
				scfKPointSplit[kpoint] = str(int(numpy.ceil(float(scfKPointSplit[kpoint])/float(factor[kpoint]))))


                        scfKPointSplit[3]=str(int(scfKPointSplit[3]))
                        scfKPointSplit[4]=str(int(scfKPointSplit[4]))
                        scfKPointSplit[5]=str(int(scfKPointSplit[5]))

			newKPointString = ' '.join(scfKPointSplit)
                        new_mod='{automatic}'
                        if scfKPointSplit[0]=="1" and scfKPointSplit[1]=="1" and scfKPointSplit[2]=="1":
                            newKPointString=""
                            new_mod="{gamma}"

                        splitInput['K_POINTS']['__content__']=newKPointString
                        splitInput['K_POINTS']['__modifier__']=new_mod
                        inputfile = AFLOWpi.retr._joinInput(splitInput)

                    return inputfile

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)



def _project_phDOS(oneCalc,ID):
    fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.phDOS.eig"%ID)
    with open(fn,"r") as fo:
        fs=fo.read()

    q_list = [k.strip() for k in re.findall("q\s*=(.*)\n",fs)]
    num_q  = len(q_list)

    re_qp=re.compile("freq.*=.*=\s*([-.0-9]+).*\n((?:\s*\(\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*\)\s*\n)+)")
    #
    qs=re_qp.findall(fs)
    by_eig=[]
    for i in qs:
        eig=[j for j in i[1].split('\n') if len(j.strip())!=0]
        qs=[map(float,k.split()[1:-1]) for k in eig]
        eig_val=[float(i[0])]
        for l in qs:
    #        x=numpy.complex(l[0],l[1])
    #        y=numpy.complex(l[2],l[3])
    #        z=numpy.complex(l[4],l[5])
            x=numpy.sqrt(l[0]**2.0+l[1]**2.0)
            y=numpy.sqrt(l[2]**2.0+l[3]**2.0)
            z=numpy.sqrt(l[4]**2.0+l[5]**2.0)
            eig_val.append([x,y,z])
        by_eig.append(eig_val)
    #by_eig=numpy.array(by_eig)
    #num_q = by_eig.shape[0]
    #num_eig = (by_eig.shape[1]-1)/3

    by_atom=[]

    for j in range(len(by_eig)):
        each_freq=[by_eig[j][0]]
        for i in range(1,len(by_eig[j])):
            each_freq.append(by_eig[j][i][0]**2+by_eig[j][i][1]**2+by_eig[j][i][2]**2)
        by_atom.append(each_freq)
    by_atom=numpy.array(by_atom)
    eig_per_q= len(by_atom)/num_q

    appdos_fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.eig.ap"%ID)


    with open(appdos_fn,"w") as fo:
        labs = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])
        ip='      q1          q2          q3            freq '+' '.join([x.rjust(15) for x in labs])+'\n'
        fo.write(ip)


        for i in range(len(by_atom)):
            ip='%32s'%q_list[i/eig_per_q]

            entries = ['%16.10f'%x for x in map(float,by_atom[i])]
            for j in entries:
                ip+=j
            ip+="\n"
            fo.write(ip)



def _sum_aproj_phDOS(appdos_fn,de=1.0):

    data,summed_weights,species = AFLOWpi.run._get_ph_weights(appdos_fn)

    data_min=numpy.min(data[:,0])*1.2
    data_max=numpy.max(data[:,0])*1.2



    num_bins=(data_max-data_min)/de
    new_bins=numpy.linspace(data_min, data_max, num=num_bins, endpoint=True)
 

    
    
    #do a new array for the bins with density,total,then a col for each species
    dos_data = numpy.zeros([len(new_bins)-1,len(species)+2])

    total_phDOS = numpy.histogram(data[:,0], bins=new_bins, )

    dos_data[:,0] = total_phDOS[1][:-1]



#    dos_data[:,1] = total_phDOS[0]

    for spec in range(1,len(species)+1):

        projected_phDOS = numpy.histogram(data[:,0], bins=new_bins, weights=summed_weights[:,spec])[0]
        dos_data[:,spec+1]=projected_phDOS
        dos_data[:,1]+=projected_phDOS

    with open(appdos_fn[:-6]+'aphdos',"w") as fo:
        fo.write('freq'.ljust(21))
        fo.write('Total'.ljust(21))

        spec = ' '.join([x.ljust(20) for x in species])
        fo.write(spec)
        for i in dos_data:
            to_write = '\n'+' '.join([x.ljust(20) for x in map(str,i.tolist())])
            fo.write(to_write)
#        for col_of_spec in col:
            
#        weights[spec]=



def _get_ph_weights(appdos_fn):

    with open(appdos_fn,"r") as fo:
        data = fo.read()

    data=data.split('\n')
    at_labels = data[0].split()[4:]

    data = numpy.asarray([map(float,x.split()[3:]) for x in data[1:] if len(x.strip())!=0])
    species=list(set(at_labels))

    species=list(collections.OrderedDict.fromkeys(at_labels))

    col_of_each_spec=collections.OrderedDict()
    for spec in species:
        for ind in range(len(at_labels)):
            if spec==at_labels[ind]:
                try:
                    col_of_each_spec[spec].append(ind+1)
                except:
                    col_of_each_spec[spec]=[ind+1]



    num_entries = data.shape[0]
    num_spec=len(species)
    summed_weights = numpy.zeros([num_entries,num_spec+1])

    summed_weights[:,0]=data[:,0]

    for spec in range(len(species)):
        cols = col_of_each_spec[species[spec]]

        
        for spec_col in cols:        
            
            summed_weights[:,spec+1]+=data[:,spec_col]

    return data,summed_weights,species
