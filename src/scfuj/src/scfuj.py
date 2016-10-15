import os
import re
import logging
import AFLOWpi
import copy 
import numpy as np
import atexit
import math 
import __main__
import itertools as it
import shutil 
import string
import itertools
import collections
import csv

def _want_txt_to_bin(fpath,fname):

    fns=fname.split(".")
    bin_file = os.path.join(fpath,fns[0]+".binary")
    if not os.path.exists( os.path.join(fpath,fname)):
        return

    fin   = open(fpath+'/'+fname,"r")
    ret=np.asarray(list(csv.reader(fin, delimiter=' ',skipinitialspace=True,
                                       quoting=csv.QUOTE_NONNUMERIC)),dtype=np.float32)
    fin.close()

    bin_data =  ret[:,0]+1j*ret[:,1]

    bout   = open(bin_file,"wb")
    np.save(bout,bin_data)
    bout.close()






def chk_species(elm):
    
    species_Nms = {#d elements
                        'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                        'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                        'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg','Sc',
                        'Ga', 'In','Y',
                       #p elements
                        'C', 'N', 'O', 'Se', 'S', 'Te','Sn','B','F','Al','Si','P','Cl',
                        'Ge','As','Br','Sb','I','Tl','Pb','Bi','Po','At','Ba',
                       #s elements
                        'H', 'Sr','Mg','Li','Be','Na','K','Ca','Rb','Cs'}

    if elm in species_Nms: return True 
    else: return False 
	

	

###################################################################################################################################################

def _getPAOfilename(atom,PAOdir=None):
        """
        Get the pseudopotential filename for the specific atomic species in input
        
        Arguments:
         - atom -- a string designating the atomic species you want the pseudofile name for
        
        Keyword Arguments:
         - pseudodir -- the path of the directory containing pseudofiles
        """
        if PAOdir==None:
            PAOdir= AFLOWpi.prep._ConfigSectionMap('prep','pao_dir')        
            if os.path.isabs(PAOdir) == False:
                configFileLocation = AFLOWpi.prep._getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                paodir =  os.path.normpath(os.path.join(configFileLocation, PAOdir))

        else:
            if os.path.isabs(PAOdir) == False:
                paodir =  os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), paodir))


        atom0=''
        atom1=''
        atom0 = atom[0]
        try:
            atom1 = atom[1]
        except:
            atom1 = ''
        species='['+atom0.upper()+atom0.lower()+']'+'['+atom1.upper()+atom1.lower()+']'
        regexFromConfig="re.compile(r'^'+species+'[._-]', re.I)"
        rex1 = eval(regexFromConfig)
        for l in os.listdir(PAOdir):
                if rex1.search(l):
			PAOfilename = l
                        logging.debug('Exiting getPAOfilename')
                        return PAOfilename
        if atom.strip()!='!':
            logging.error('Missing PAO File')
            logging.debug('Exiting getPAOfilename')
        return None

###################################################################################################################################################



#def scfprep(allAFLOWpiVars,refFile,pseudodir=None,paodir=None,build_type='product'):
def scfprep(calcs,paodir=None):
    output_calcs=copy.deepcopy(calcs)
    for ID,oneCalc in calcs.iteritems():
        AFLOWpi.prep._addToBlock(oneCalc,ID,'PREPROCESSING','oneCalc,ID = AFLOWpi.scfuj._oneScfprep(oneCalc,ID)') 
    return output_calcs

def _oneScfprep(oneCalc,ID,paodir=None):
	"""
	Read a ref file (str), the dictionary defining the calculations, and the dictionary with the aflowkeys Create the dir 
        tree for the run and return a dictionary dictAllcalcs with all the calculations. Store the dictionary in a log file

        Arguments:
         - allAFLOWpiVars -- all the variables that you want to make a list of combinations of calculations from
         - refFile    -- string that contains the input ref file file path
        Keyword Arguments:
         - pseudodir  -- path of the directory that contains your Pseudopotential files
         - paodir     -- path of the directory that contains pseudo atomic orbital files for the respective pseudo-potentials
	 - calcs      -- Dictionary of dictionaries of calculations - this can be empty if it is the initial acbn0 run.
	
	"""

        if 'SCFUJ Iteration' in  oneCalc['__status__'].keys():
            return oneCalc,ID

        configFileLocation = AFLOWpi.prep._getConfigFile()
        if not os.path.exists(configFileLocation):
            configFileLocation = oneCalc['_AFLOWPI_CONFIG_']
#            break
        
	output_oneCalc = copy.deepcopy(oneCalc)
	if paodir==None:
		paodir = AFLOWpi.prep._ConfigSectionMap('prep','pao_dir')
                if os.path.isabs(paodir) == False:
                    configFileLocation = os.path.dirname(configFileLocation)
                    paodir =  os.path.join(configFileLocation, paodir)

        else:
            if os.path.isabs(paodir) == False:
                paodir =  os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), paodir))
	
        oneCalc = copy.deepcopy(oneCalc)



        #Modify inputfile to include Hubbard parameters
        inputfile = oneCalc['_AFLOWPI_INPUT_']
        #check if u vals are already in input and start with those.



        #Assign initial U values
#        species=list(set(AFLOWpi.retr._getPosLabels(inputfile)))

        species = re.findall("(\w+).*UPF",inputfile)
        splitInput = AFLOWpi.retr._splitInput(inputfile)

        Uvals = {}
        uppered_system={}
        for k,v in splitInput['&system'].iteritems():
            uppered_system[k.upper()]=v

        for isp in range(len(species)):
            if 'HUBBARD_U(%s)' % (isp+1) in uppered_system.keys():
                startingUval = float(uppered_system['HUBBARD_U(%s)' % (isp+1)])

                Uvals[species[isp]] = startingUval

            else:
                Uvals[species[isp]] = 0.001	


	AFLOWpi.prep._modifyVarVal(oneCalc,ID,varName='uValue',value=Uvals)
        new_calc = updateUvals(oneCalc,Uvals,ID=ID)

        new_inputfile = new_calc['_AFLOWPI_INPUT_']
        calc_label = ID

        output_onecalc = new_calc
        output_oneCalc['_AFLOWPI_INPUT_'] = new_inputfile
        try:
            output_oneCalc['__status__']['SCFUJ Iteration']=0
        except:
            pass

	

	#adding config file location to the calculations
	configFile = AFLOWpi.prep._getConfigFile()

        maketree(output_oneCalc,calc_label, paodir=paodir)


        """Generate a new uValLog.log file to wipe away if there is an old one"""
        
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_uValLog.log' % ID),'w') as uValLogFile:
            uValLogFile.write('')

        
	return output_oneCalc,calc_label

def updateUvals(oneCalc, Uvals,ID=None):
	"""
	Modify scf input file to do a lda+u calculation.
	
	Arguments:
	 - oneCalc      -- Dictionary of one calculation 
	 - Uvals	-- Dictionary of Uvals

	"""
	try:
		inputfile = oneCalc['_AFLOWPI_INPUT_']
		#Get species
#                species=list(set(AFLOWpi.retr._getPosLabels(inputfile)))
		species = re.findall("(\w+).*UPF",inputfile)

                inputDict = AFLOWpi.retr._splitInput(inputfile)
                inputDict['&system']['lda_plus_u']='.TRUE.'
                for isp in range(len(species)):
                    hub_entry = 'Hubbard_U(%s)'% str(isp+1)
                    inputDict['&system'][hub_entry] = Uvals[species[isp]]
         	#Update inputfile
                oneCalc['_AFLOWPI_INPUT_']=AFLOWpi.retr._joinInput(inputDict)

	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)


	return oneCalc


###################################################################################################################################################

def maketree(oneCalc,ID, paodir=None):

	"""
        Make the directoy tree and place in the input file there
        
        Arguments:
         - calcs -- Dictionary of dictionaries of calculations

        Keyword Arguments:
         - pseudodir   -- path of pseudopotential files directory 
	 - paodir	- path of pseudoatomic orbital basis set
        """

#        for ID,oneCalc in calcs.iteritems():
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % ID),'w') as scfujInfileObj:
            scfujInfileObj.write(oneCalc['_AFLOWPI_INPUT_'])

        '''save the calc just in case so it's updated witht he scfuj U value info in the input file as well as _AFLOWPI_INPUT_'''

#        if AFLOWpi.prep._findInBlock(oneCalc,ID,'LOADCALC','''oneCalc = AFLOWpi.prep._loadOneCalc('%s','%s')''' % (oneCalc['_AFLOWPI_FOLDER_'],ID) )==False:

#            AFLOWpi.prep._addToBlock(oneCalc,ID,'LOADCALC','''try:
#   oneCalc = AFLOWpi.prep._loadOneCalc('%s','%s')
#except:
#   pass
#''' % (oneCalc['_AFLOWPI_FOLDER_'],ID) )            


	#Get wanT directory from config file and copy bands.x from it to tre
	wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
        if os.path.isabs(wantdir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            wantdir =  os.path.join(configFileLocation, wantdir)

	wantBandsExec = os.path.join(wantdir,'bands.x')
        if os.path.exists(wantBandsExec)==False:
            print 'WanT bands executable not found. Check your config file to make sure want_dir path is correct and that bands.x is in that directory..Exiting'
            logging.error('WanT bands executable not found. Check your config file to make sure want_dir path is correct and that bands.x is in that directory..Exiting')
            raise SystemExit


	wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
        if os.path.isabs(wantdir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            wantdir =  os.path.join(configFileLocation, wantdir)

	espressodir = AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if os.path.isabs(espressodir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            espressodir =  os.path.join(configFileLocation, espressodir)

	pdosExec = 'projwfc.x'
	pdosExec = os.path.join(espressodir,pdosExec)
        if os.path.exists(pdosExec)==False:
            print 'ProjectWFC executable not found. Check your config file to make sure engine_dir path is correct and that projwfc.x is in that directory..Exiting'
            logging.error('ProjectWFC executable not found. Check your config file to make sure engine_dir path is correct and that projwfc.x is in that directory..Exiting')
            raise SystemExit
	dosExec = 'dos.x'
	dosExec = os.path.join(espressodir,dosExec)
        if os.path.exists(dosExec)==False:
            print 'DOS executable not found. Check your config file to make sure engine_dir path is correct and that dos.x is in that directory..Exiting'
            logging.error('DOS executable not found. Check your config file to make sure engine_dir path is correct and that dos.x is in that directory..Exiting')
            raise SystemExit
	scfExec = 'pw.x'
	scfExec = os.path.join(espressodir,scfExec)
        if os.path.exists(scfExec)==False:
            print 'SCF executable not found. Check your config file to make sure engine_dir path is correct and that pw.x is in that directory..Exiting'
            logging.error('SCF executable not found. Check your config file to make sure engine_dir path is correct and that pw.x is in that directory..Exiting')
            raise SystemExit
            


	try:
		if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
			AFLOWpi.prep.totree(wantBandsExec,{ID:oneCalc},rename='want_bands.x')
			AFLOWpi.prep.totree(pdosExec,{ID:oneCalc})
			AFLOWpi.prep.totree(dosExec,{ID:oneCalc})
			AFLOWpi.prep.totree(scfExec,{ID:oneCalc})
                else:
			AFLOWpi.prep.totree(wantBandsExec,{ID:oneCalc},rename='want_bands.x',symlink=True)
			AFLOWpi.prep.totree(pdosExec,{ID:oneCalc},symlink=True)
			AFLOWpi.prep.totree(dosExec,{ID:oneCalc},symlink=True)
			AFLOWpi.prep.totree(scfExec,{ID:oneCalc},symlink=True)
	except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
##############################################################################################################
	#move acbn0.py to dir tree
	acbn0Path = os.path.join(AFLOWpi.__path__[0],'scfuj','acbn0_support', 'acbn0.py')
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower() == 'false':
		AFLOWpi.prep.totree(acbn0Path,{ID:oneCalc},symlink=True)
	else:
		AFLOWpi.prep.totree(acbn0Path,{ID:oneCalc},symlink=False)

	#move integs.pyc to dir tree
	integPath = os.path.join(AFLOWpi.__path__[0],'scfuj','acbn0_support', 'integs.py')	
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower() == 'false':
		AFLOWpi.prep.totree(integPath,{ID:oneCalc},symlink=True)
	else:
		AFLOWpi.prep.totree(integPath,{ID:oneCalc})

	#move integs.pyc to dir tree
	pyintsPath = os.path.join(AFLOWpi.__path__[0],'scfuj','acbn0_support', 'pyints.py')	
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower() == 'false':
		AFLOWpi.prep.totree(pyintsPath,{ID:oneCalc},symlink=True)
	else:
		AFLOWpi.prep.totree(pyintsPath,{ID:oneCalc})



        moleculePath = os.path.join(AFLOWpi.__path__[0],'scfuj','acbn0_support', 'Molecule.py')
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower() == 'false':
            AFLOWpi.prep.totree(moleculePath,{ID:oneCalc},symlink=True)

	else:

		AFLOWpi.prep.totree(moleculePath,{ID:oneCalc})

############################################################################################################
	#Get paodir from config file
        try:
            if paodir==None:
		paodir= AFLOWpi.prep._ConfigSectionMap('prep','pao_dir')
                if os.path.isabs(paodir) == False:
                    configFileLocation = AFLOWpi.prep._getConfigFile()
                    configFileLocation = os.path.dirname(configFileLocation)
                    paodir =  os.path.join(configFileLocation, paodir)

            else:

                if os.path.isabs(paodir) == False:
                    paodir =  os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), paodir))

                else:
                    pass
#                    print 'Can not find PAO file. Exiting'
#                    logging.error('Can not find PAO file. Exiting')
#                    raise SystemExit
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            print 'Can not find PAO file. Exiting'
            logging.error('Can not find PAO file. Exiting')
            raise SystemExit




        #Copy PAO files 

        try:
            work_dir = oneCalc['_AFLOWPI_FOLDER_'] 
            for key in oneCalc.keys():
                    v = oneCalc[key]
                    if re.search(r'_AFLOWPI_[A-Z][0-9]*_', key):
                            vp = AFLOWpi.scfuj._getPAOfilename(v.strip('0123456789'),paodir)		


                            try:
                                    a = os.path.join(paodir,vp)
                            except AttributeError:
                                if v=='!':
                                    continue
                                logging.debug('Cannot find correct PAO files in %s ...Exiting' % paodir)
                                print 'cannot find correct PAO files in %s ...Exiting' % paodir
                                raise SystemExit

                            newPAOFileNm = os.path.join(work_dir,v.strip('0123456789')+"_basis.py")
                            print 'Copying '+a+' to '+ newPAOFileNm+'\n'
                            logging.info('Copying '+a+' to '+ newPAOFileNm)
                            if AFLOWpi.prep._ConfigSectionMap('prep','copy_pseudos').lower() == 'false':
                                    try:
                                            os.symlink(a,newPAOFileNm)
                                    except OSError:
                                        os.system('rm -fr %s' % newPAOFileNm)
                                        try:
                                            os.symlink(a,newPAOFileNm)
                                        except:
                                            logging.error('cant copy PAO!')
                                            raise SystemExit
                            else:
#                                    print a
#                                    print newPAOFileNm

                                    if not os.path.exists(newPAOFileNm):
                                        try:
                                            shutil.copy(a, newPAOFileNm)
                                        except:
                                            try:
                                                os.system('rm -fr %s' % newPAOFileNm)
                                                shutil.copy(a, newPAOFileNm)
                                            except:
                                                pass
                                    else:
                                        pass

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


        '''save the inital list of things to exec so first iteration it will run through'''

#        oneCalc['__runList__']=[]
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        return oneCalc

def nscf_nosym_noinv(oneCalc,ID=None,kpFactor=1.50,unoccupied_states=False):
        """
        Add the ncsf input to each subdir and update the master dictionary
        
        Arguments:
         - calc_copy -- dictionary of one calculation
         - kpFactor -- multiplicative factor for kpoints from SCF to DOS (default: 2)
        """

	calc_copy = copy.deepcopy(oneCalc)
        try:
                subdir = calc_copy['_AFLOWPI_FOLDER_']
		a = ID+'.in'
                inputfile = oneCalc['_AFLOWPI_INPUT_']
                '''we need nscf and not bands because we need output for HOMO'''
#                inputDict=AFLOWpi.retr._splitInput(inputfile)

#                inputfile=AFLOWpi.retr._joinInput(inputDict)
                a = ID+'.out'
                splitInput = AFLOWpi.retr._splitInput(inputfile)

                try:

                        for ID_old in oneCalc['prev']:
                            try:
                                a = ID_old+'.out'
                                outfile = open(os.path.join(subdir,a),'r').read()
                                match1 = re.search(r'Kohn-Sham states=',outfile)
                                if match1==None:
                                    continue
                                m1 = match1.end()
                                m2 = m1 + 15
                                nbnd = int(1.75*int(outfile[m1:m2]))
                                if unoccupied_states==True:
                                    splitInput['&system']['nbnd']= nbnd                        
                                else:
                                    try:
                                        del splitInput['&system']['nbnd']
                                    except:
                                        pass
                                print 'Number of bands to be Calculated %s: '% nbnd
                                logging.info('Number of bands to be Calculated %s: '% nbnd)
                                break
                            except Exception,e:
                                pass

			'''checks to see if nbnd has already been added to the file - in any case add/replace nbnd'''
#                        if "occupations" in splitInput['&system'].keys():
#                            if splitInput['&system']["occupations"]!= 'smearing':
#                                splitInput['&system']["occupations"]='"tetrahedra"'
#                        else:
#                            splitInput['&system']["occupations"]='"tetrahedra"'

			'''Add nosym = .true. and noinv = .true. to file '''

                        try:
                            splitInput['&electrons']['conv_thr']='1.0D-8'
                            splitInput['&system']['nosym']='.true.'
                            splitInput['&system']['noinv']='.true.'

                            splitInput['&control']['calculation']='"nscf"'
                            inputfile=AFLOWpi.retr._joinInput(splitInput)
                            '''writes an input for band_plot.x to process the correct number of bands calculated'''
                        
                            
                        except Exception,e:
				AFLOWpi.run._fancy_error_log(e)
                                print e

                except Exception,e:
                        AFLOWpi.run._fancy_error_log(e)
                        print 'SCF outputfile not found: nbnd is default'
                        logging.info('SCF outputfile not found: nbnd is default')
                        pass
                '''checks to see if "crystal_b" has already been substituted for "automatic"'''

                try:

                    mod = splitInput['K_POINTS']['__modifier__'].upper()
                    mod=mod.strip('{}()')
                    if mod=='GAMMA':
                        splitInput['K_POINTS']['__content__']='2 2 2 0 0 0'
                        splitInput['K_POINTS']['__modifier__']='{automatic}'
                        inputfile = AFLOWpi.retr._joinInput(splitInput)
                    else:
                        splitInput=AFLOWpi.retr._splitInput(inputfile)
                        scfKPointString = splitInput['K_POINTS']['__content__']
			scfKPointSplit = [float(x) for x in scfKPointString.split()]

                        before_kpf=reduce(lambda x, y: x*y, scfKPointSplit )
                        scaling_kpf=before_kpf/150.0
                        if scaling_kpf>1.0:
                            kpFactor=1

                            
                            

			for kpoint in range(len(scfKPointSplit)-3):
				scfKPointSplit[kpoint] = str(int(math.ceil(scfKPointSplit[kpoint]*kpFactor)))
			for kpoint in range(3,len(scfKPointSplit)):
				scfKPointSplit[kpoint] = '0'

			newKPointString = ' '.join(scfKPointSplit)

                        splitInput['K_POINTS']['__content__']=newKPointString
                        splitInput['K_POINTS']['__modifier__']='{automatic}'
                        inputfile = AFLOWpi.retr._joinInput(splitInput)

                except Exception,e:
			AFLOWpi.run._fancy_error_log(e)


		calc_label = ID + "_nscf"

                a = calc_label+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile)
                new_inputfile.close()

                output_calc = calc_copy
                output_calc['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]
        except IOError as e:
                logging.error("%s not in %s" % (ID,calc_copy['_AFLOWPI_FOLDER_']))

                
        return output_calc,calc_label


def projwfc(oneCalc,ID=None):
	'''
	Run projwfc on each calculation

	Arguments:
	 - oneCalc -- dictionary of a single calculation

	'''

        temp_dir= AFLOWpi.prep._get_tempdir()
	calc_copy=copy.deepcopy(oneCalc)
        output_calc = {}
        prefix = oneCalc['_AFLOWPI_PREFIX_']
        try:
                subdir = oneCalc['_AFLOWPI_FOLDER_']
        
                try:
                    prefix = AFLOWpi.retr._prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])


                except Exception,e:

                    prefix = oneCalc['_AFLOWPI_PREFIX_']
                

		inputfile = """&PROJWFC
  prefix='%s'
  filpdos='./%s_acbn0'
  outdir='%s'
  lwrite_overlaps=.TRUE.
  lbinary_data  = .TRUE.
/
"""%(prefix,ID,temp_dir)

		calc_label = ID + '_pdos'		

                a = calc_label+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile)
                new_inputfile.close()

                output_calc[calc_label] = calc_copy
                output_calc[calc_label]['_AFLOWPI_INPUT_'] = inputfile


        except IOError as e:
                logging.error("%s not in %s" % (ID,oneCalc['_AFLOWPI_FOLDER_']))
	
	return output_calc[calc_label], calc_label


###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

def WanT_bands(oneCalc,ID=None,eShift=15.0,num_points=1000,cond_bands=0):
	'''
		Make input files for  WanT bands calculation

		Arguments:
        	 - calc_copy -- dictionary of dictionaries of calculations
	'''
	output_calc = {}
        temp_dir= AFLOWpi.prep._get_tempdir()
	calc_copy=copy.deepcopy(oneCalc)
	subdir = calc_copy['_AFLOWPI_FOLDER_']
	nspin = int(AFLOWpi.scfuj.chkSpinCalc(calc_copy,ID))

	### GET KPATH for respective symmetry #####
	special_points, band_path = AFLOWpi.retr._getHighSymPoints(calc_copy)
	lblList = []; kpathStr = ""
	for i in band_path.split('|'):
                lblList += i.split('-')
        for i in lblList:
                kpathStr += "%s\t%.5f\t%.5f\t%.5f\n"%(i,special_points[i][0],special_points[i][1],special_points[i][2])
	############################################ 
        try:
            nscf_ID=ID+'_nscf'
            Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
            eShift=float(Efermi)+10.0
        except:
            eShift=10.0

        eShift=5.0
            
        try:
            prefix = AFLOWpi.retr._prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        except:
            prefix = oneCalc['_AFLOWPI_PREFIX_']
                
#        if temp_dir!='./':
#            subdir=temp_dir

        nbnd_string= ''
        if cond_bands!=0 and type(cond_bands) == type(452):
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False))+cond_bands
                nbnd_string = 'atmproj_nbnd    = %s' % nbnds
	elif cond_bands==True: 
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=True))
                nbnd_string = 'atmproj_nbnd    = %s' % nbnds
	elif cond_bands==0 or cond_bands==False:
#            nbnd_string=""
            nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False))
            nbnd_string = 'atmproj_nbnd    = %s' % nbnds
        
#        if nbnd==None:

 #           else:
#                nbnd=int(nbnd)
#                nbnd_string = 'atmproj_nbnd    = %s' % nbnd


	inputfile = """ &INPUT							  
prefix 		= '%s_TB' 		  	                        
postfix 	= \'_WanT\'		        
work_dir	= \'./\'		 
datafile_dft	= \'./%s_TB.save/atomic_proj.dat\'	 
nkpts_in        = %d                 
nkpts_max	= %s		 
do_orthoovp	= .FALSE.
atmproj_sh      = %f                         
%s
"""%(ID,ID,len(lblList),num_points,eShift,nbnd_string)
#"""%(prefix,subdir,subdir,prefix,len(special_points),eShift,nbnd_string)

	if nspin == 1:
		inputfile += "fileout	= '%s_bands_want.dat'\n/\n"%ID
		# Insert k-path 
		inputfile += kpathStr

		calc_label = ID + "_WanT_bands"
		
	        a = calc_label+'.in'
	        new_inputfile = open(os.path.join('./',a),'w')
	        new_inputfile.write(inputfile)
	        new_inputfile.close()

	        output_calc = calc_copy
	        output_calc['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]

		single_output_calc = {calc_label:output_calc}

		return single_output_calc

	else:
		inputfile1 = inputfile
		#Create two input files for spin-up and spin-down components.
		inputfile +=  """fileout	= '%s_bands_want_up.dat'  \nspin_component	  =    "up"\n/\n"""%ID
		inputfile1 += """fileout	= '%s_bands_want_down.dat'\nspin_component   =    "down"\n/\n"""%ID
		#Insert k-path
		inputfile += kpathStr
		inputfile1 += kpathStr

                calc_label_up = ID + "_WanT_bands_up" 
                a = calc_label_up+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile)
                new_inputfile.close()
                output_calc_up= calc_copy
                output_calc_up['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]


                calc_label_down = ID + "_WanT_bands_down" 
                a = calc_label_down+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile1)
	        new_inputfile.close()
                output_calc_down = calc_copy
                output_calc_down['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]


		output_calc = {calc_label_up:output_calc_up, calc_label_down:output_calc_down}		
		return output_calc


###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

def chkSpinCalc(oneCalc,ID=None):

	'''
		Check whether an calculation is spin polarized or not.
	
		Arguments:
		
		--oneCalc : dictionary of a single calculation.

	'''
        try:
            nspin = int(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['nspin'])
            return nspin
        except Exception,e:
            return 1

        oneCalcID = ID#oneCalc['_AFLOWPI_PREFIX_'][1:]
        subdir = oneCalc['_AFLOWPI_FOLDER_']
        scfOutput = '%s.out' % oneCalcID

        regex = re.compile(r"(spin.*)\n",re.MULTILINE)


        for ID_old in oneCalc['prev']:
            try:
                fin = os.path.join(subdir,ID_old+'.out')
                with open(fin,'r') as fin_obj:
                    lines = fin_obj.read()

                if len(regex.findall(lines)) != 0:
                    return 2
                else:
                    return 1

            except Exception,e:
                pass


	

def evCurveMinimize(calcs,config=None,pThresh=10.0,final_minimization = 'vc-relax'):

	for ID,oneCalc in calcs.iteritems():
		execString = '''
newOneCalc = AFLOWpi.scfuj. _oneMinimizeCalcs(oneCalc, ID, pThresh=%f)
AFLOWpi.prep._saveOneCalc(oneCalc,ID)
	'''%(pThresh)
		AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN', execString)
	AFLOWpi.run._skeletonRun(calcs)
	#if final_minimization != None:
#		calcs = AFLOWpi.prep.changeCalcs(calcs, 'calculation', final_minimization)
#		AFLOWpi.run.scf(calcs)
#		calcs = AFLOWpi.prep.changeCalcs(calcs, 'calculation', 'scf')
#	AFLOWpi.run.scf(calcs)

	return calcs

def _oneMinimizeCalcs(oneCalc,ID,config=None,pThresh=10.0):
	'''
		Get equilibrium volume using evfit to Murnaghan's EOS with 5 volumes in +/-20% of input volume
	'''



	__submitNodeName__ =  __main__. __submitNodeName__
	execPrefix = ''
	execPostfix = ''
	if config!=None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug('forced config %s' % config)
	else:
		try:
			config = AFLOWpi.prep._getConfigFile()
			AFLOWpi.prep._forceGlobalConfigFile(config)
		except Exception,e:
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
	

	def runCalcList(alatList,oneCalc,newCalc):

		S="Begin running calculations for lat para " + str(alatList) + "of %s"%oneCalc['_AFLOWPI_FOLDER_']
		logging.info(S)
		engList = []; stressList = []
		for i in alatList:
			try:
				subdir = oneCalc['_AFLOWPI_FOLDER_']
				inputfile = oneCalc['_AFLOWPI_INPUT_']

				inputfile = celldmRE.sub("celldm(1)=%f"%i,inputfile)
				newCalc['_AFLOWPI_INPUT_']= inputfile

				a = ID + '.in' 
				new_inputfile = open(os.path.join(subdir,a),'w')
				new_inputfile.write(inputfile)
				new_inputfile.close()
			
					
				AFLOWpi.run._oneRun(__submitNodeName__,newCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)		
				
				outFile = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)
				outFileString = file(outFile,'r').read()
				energyRegex = re.compile(r'(?:(?:(?:(?:\!\s+)total)|(?:Final)) en\w+\s*=\s+(.+?)Ry)',re.MULTILINE)
				finalEnergy=energyRegex.findall(outFileString)
				stressRegex = re.compile(r'P=\s*(.*)\n')
				finalStress = stressRegex.findall(outFileString)
				

				if len(finalEnergy):
					energy  = float(finalEnergy[-1]); engList.append(energy)
					stress = float(finalStress[0]);	stressList.append(stress)
					logging.info("E-V curve point: %f a.u. %f Ry" %(i, energy))
					print "E-V curve point: %f a.u. %f Ry" %(i, energy)	
				else:
					logging.error("ERROR: E-V curve at celldm(1) = %f FAILED"%i)
					raise SystemExit

			except Exception, e:
	                        AFLOWpi.run._fancy_error_log(e)

		return zip(stressList, engList, alatList)

	def turningPts(tmplst):
                lstarr = np.array(tmplst)
                de = np.diff(lstarr.T[1])
                npts = np.sum(de[1:] * de[:-1] < 0)
                logging.info("%d turning points found in E-V curve"%npts)
                return npts


        subdir = oneCalc['_AFLOWPI_FOLDER_']
	oneCalc['_AFLOWPI_CONFIG_']=config
	
        newCalc = oneCalc
	inputfile = newCalc['_AFLOWPI_INPUT_']


	ibravRE = re.compile(r"ibrav.*?(\d+)")
        ibrav = int(float(ibravRE.findall(inputfile)[0])) 	

	

	if ibrav == 1 or ibrav == 2 or ibrav ==3:
		logging.info("Entering structural optimization via 5 point E-V curve")
		try:
			celldmRE = re.compile('celldm\(1\).*?=.*?([0-9]*\.?[0-9]+),?')	
			celldm1 = float(celldmRE.findall(inputfile)[0])
			#Get list of lattice parameters in the range 70% V0 to 130% V0
			alatList = list(np.arange((celldm1**3*.7)**(1./3.),(celldm1**3*1.3)**(1./3.), ((celldm1**3*1.2)**(1./3.)-(celldm1**3*.8)**(1./3.))/5))
			#Get first set of datapoints for E-V curve
			dataList = runCalcList(alatList,oneCalc,newCalc)
			
                        #Check whether a turning point is reached in E-V curve
                        tmpList = dataList

                        nTurnPts = turningPts(tmpList) 

                        #Sort on stress
                        tmpList.sort(key=lambda tup: tup[0])
                        #tmpfout = file("tmpEV.txt",'a')
			while nTurnPts < 1:
                                
                         #       for i in tmpList:
                         #               tmpfout.write("%f\t %f\t %f\n"%(i[2], i[1], i[0]))
                         #       tmpfout.close()

				#Check if volume expansion or compression is needed to arrive at V0
				if tmpList[0][0] > pThresh:
					celldm1 = dataList[0][2]*1.1
				else:
					celldm1 = dataList[0][2]*0.9

				#Get list of lattice parameters in the range 70% V0 to 130% V0
				alatList = list(np.arange((celldm1**3*.7)**(1./3.),(celldm1**3*1.3)**(1./3.), ((celldm1**3*1.2)**(1./3.)-(celldm1**3*.8)**(1./3.))/5))
				#Get new set of points
				newdataList = runCalcList(alatList,oneCalc,newCalc)
				#Add to old set of points to make curve better
				dataList = dataList + newdataList
                                #Find the number of turning points
                                nTurnPts = turningPts(dataList)

				tmpList = newdataList                               
                                tmpList.sort(key=lambda tup: tup[0])
				

                        evinFile = os.path.join(subdir,'%s_evx.in'%ID)
                        evoutFile = os.path.join(subdir,'%s_evx.out'%ID)
                        fout = file(evinFile,'a')
			dataList.sort(key=lambda tup: tup[2])
                        #Write data to file
                        for i in dataList:
                                fout.write("%f\t %f\n"%(i[2], i[1]))
                        fout.close()

			engineDir=AFLOWpi.prep._ConfigSectionMap("prep","engine_dir")
			evfitPath=os.path.join(engineDir,'ev.x')
			evfitString = "cat<<! | %s \nau\nsc\n4\n%s\n%s\n!\n"%(evfitPath,evinFile,evoutFile)
			logging.info("Starting fitting E-V curve with Murnaghan's EOS")
			print "Starting Fitting E-V curve with Murnaghan's EOS"
			os.system(evfitString)
			evoutFileString = file(evoutFile,'r').read()
			logging.info("Finished fitting E-V curve with Murnaghan's EOS")
			print "Finished Fitting E-V curve with Murnaghan's EOS"
			
			a0regex = re.compile('a0.*?=.*?([0-9]*\.?[0-9]+)\s*a\.u\.')
			a0 = float(a0regex.findall(evoutFileString)[0])
			inputfile = celldmRE.sub("celldm(1)=%f"%a0,inputfile)
			newCalc['_AFLOWPI_INPUT_']= inputfile
			
		        if AFLOWpi.prep._findInBlock(oneCalc,ID,'LOADCALC','''oneCalc = AFLOWpi.prep._loadOneCalc('%s','%s')''' % (oneCalc['_AFLOWPI_FOLDER_'],ID))==False:
		            AFLOWpi.prep._addToBlock(oneCalc,ID,'LOADCALC','''oneCalc = AFLOWpi.prep._loadOneCalc('%s','%s')''' % (oneCalc['_AFLOWPI_FOLDER_'],ID) )
		        AFLOWpi.prep._saveOneCalc(newCalc,ID)
 
			#Save new inputfile
                        a = ID + '.in'
			new_inputfile = open(os.path.join(subdir,a),'w')
			new_inputfile.write(inputfile)
			new_inputfile.close()


		except Exception, e:
			AFLOWpi.run._fancy_error_log(e)

		return newCalc
	else:
		print "Minimization for ibrav = %d not implemented" %ibrav
		logging.error("Minimization for ibrav = %d not implemented"%ibrav)
		raise SystemExit
def acbn0(oneCalc,projCalcID,byAtom=False):
        '''
        



        '''

	oneCalcID = '_'.join(projCalcID.split('_')[:-1])#oneCalc['_AFLOWPI_PREFIX_'][1:]	
        print oneCalcID
	subdir = oneCalc['_AFLOWPI_FOLDER_']	
	nspin = chkSpinCalc(oneCalc,oneCalcID)

	def get_orbital(elm):

                #d elements
                trM_Nms ={ 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                           'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                           'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg','Sc',
                           'Ga', 'In','Y'}

                #p elements
                pElm_Nms = {'C', 'N', 'O', 'Se', 'S', 'Te','Sn','B','F','Al','Si','P','Cl',
                            'Ge','As','Br','Sb','I','Tl','Pb','Bi','Po','At', 'Ba'}

                #s elements
                sElm_Nms = {'H', 'Sr','Mg','Li','Be','Na','K','Ca','Rb','Cs'}


	
		if elm in trM_Nms: return 2
		elif elm in pElm_Nms: return 1
		elif elm in sElm_Nms: return 0
		
	
	def gen_input(oneCalcID,subdir,nspin):
		try:
			#Get cell parameters, arranged as a single string in pattern a1i, a1j, a1k, a2i, a2j, a2k, a3i, a3j, a3k
			a,cell=AFLOWpi.retr._getCellParams(oneCalc,oneCalcID)
                        """THINK OF A BETTER WAY TO CHECK THIS"""
                        try:
                            if cell.getA()[0][0]<2.0 and cell.getA()[0][0]>0.0:
                                cellParaMatrix=a*cell

                            else:
                                cellParaMatrix=1.0*cell
                        except Exception,e:
                            print e
                        l=cellParaMatrix.tolist()
			cellParaStr = ""
                        """THINK OF A BETTER WAY TO CHECK THIS"""
			for i in range(3):
				cellParaStr += str(l[i]).strip('[]')
				if i != 2:cellParaStr += ' ,'

			#Get atomic positions in cartesian coordinates, in a single string in pattern x1,y1,z1, x2, y2, z2, ..., xn, yn, zn
			scfOutput = '%s.out' % oneCalcID
			fin = file(os.path.join(subdir,scfOutput),'r')
			lines = fin.read()

			atmPosRegex = re.compile(r"positions \(alat units\)\n((?:.*\w*\s*tau\(.*\)\s=.*\(.*\)\n)+)",re.MULTILINE) 
                        try:
                            positions = AFLOWpi.retr.getPositionsFromOutput(oneCalc,oneCalcID)
                            positions = AFLOWpi.retr._cellMatrixToString(positions)
                            
                            if positions=="":
                                splitInput = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
                                positions= splitInput["ATOMIC_POSITIONS"]["__content__"]
                        except:
                            splitInput = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
                            positions= splitInput["ATOMIC_POSITIONS"]["__content__"]


			lines1 = atmPosRegex.findall(lines)[0]
			atmPosRegex1 = re.compile(r".*=.*\((.*)\)\n+",re.MULTILINE)

			atmPos = atmPosRegex1.findall(lines1)

                        try:
                            positions=np.asarray([[float(i.split()[1]),float(i.split()[2]),float(i.split()[3])] for i in positions.split("\n") if len(i.strip())!=0])
                        except:
                            positions=np.asarray([[float(i.split()[0]),float(i.split()[1]),float(i.split()[2])] for i in positions.split("\n") if len(i.strip())!=0])
                        positions=AFLOWpi.retr._convertCartesian(positions,cellParaMatrix,scaleFactor=1.0)
                        atmPos=AFLOWpi.retr._cellMatrixToString(positions).split("\n")
                        
			atmPosList = []

			for i in atmPos:atmPosList.append(map(float,i.split()))
			atmPosStr = ""

		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)
                        raise SystemExit


                for i in range(len(atmPosList)):
			try:
		#Convert fractional atomic coordinates to cartesian coordinates using the lattice vectors (cell parameters)
#				atmPosStr += str(list(np.array(atmPosList[i])*a)).strip('[]')
				atmPosStr += str(list(np.array(atmPosList[i]))).strip('[]')
				if i != len(atmPosList)-1:
					atmPosStr += ' ,'
			except Exception,e:
				AFLOWpi.run._fancy_error_log(e)       
			try:
				#Get the list of atom labels arranged as a single string
				atmLblRegex = re.compile(r"\d+\s+(\w+).*=.*\n+",re.MULTILINE)
				atmLbls = atmLblRegex.findall(lines1)

				#Strip digits from atom lables
				for i in range(len(atmLbls)):atmLbls[i]=atmLbls[i].strip('0123456789')

				atmLblsStr= str(atmLbls).strip('[]').replace("'","")
				natoms = len(atmLblsStr.split())

				#Get the list of atom species
				atmSpRegex = re.compile(r"atomic species   valence    mass     pseudopotential\n(?:.*\(.*\)\n)+",re.MULTILINE)
				atmSpList = []
                                atmSpStrList = atmSpRegex.findall(lines)[0].split('\n')[1:-1]

                                for i in atmSpStrList:atmSpList.append(i.split()[0])
				fin.close()

				#Get list of orbitals
				projOut = projCalcID + '.out'
				fin = file(os.path.join(subdir,projOut), 'r')
				proj_lines = fin.read()

				inFileList = []

			except Exception,e:
				AFLOWpi.run._fancy_error_log(e)


                #in the case we want the eff_U from each atom we make a list of the site numbers
                #and use a slightly different regex to get the l numbers
                if byAtom:
                    atmSpList= [str(x) for x in range(1,len(atmPosList)+1)]

		#For each atomic species
		for atmSp in atmSpList:
			logging.info("Creating acbn0 inpufile for %s"%atmSp)
			try:
				#Get orbital type to apply Hubbard correction
                                try:
                                    atmSp=oneCalc['__scfuj_label_mapping__'][atmSp]
                                except:
                                    pass
                                ql = get_orbital(atmSp.strip('0123456789'))
#                                print 'ql',ql
                                if byAtom==False:

       				#Get list of all orbitals of type ql of the same species
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*\(%s.*\).*\(l=%d.*\)\n"%(atmSp.strip('0123456789'),ql),re.MULTILINE)
#                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*\(%s.*\).*\(l=%d.*\)\n"%(atmSp,ql),re.MULTILINE)
                                    
                                    eqOrbList = map(int, map(float, eqOrbRegex.findall(proj_lines)))
                                    red_basis = [x - 1 for x in eqOrbList]

				#Get ones relevant for hubbard center
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*%s.*l=%d.*\)\n"%(atmSp,ql),re.MULTILINE)
                                else:

                                    getSpecByAtomRegex = re.compile(r"state #\s*(?:\d*): atom\s*%s\s*\(([a-zA-Z]+)"%atmSp)
                                    speciesFromNum = getSpecByAtomRegex.findall(proj_lines)[-1]
                                    ql = get_orbital(speciesFromNum.strip('0123456789'))

                                    #Get list of all orbitals of type ql of the same atom
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*%s.*l=%d.*\n"%(speciesFromNum,ql),re.MULTILINE)          
#                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*%s.*l=%d.*\n"%(speciesFromNum.strip('0123456789'),ql),re.MULTILINE)          


                                    eqOrbList = map(int, map(float,eqOrbRegex.findall(proj_lines)))
                                    red_basis = [x - 1 for x in eqOrbList]

				#Get ones relevant for hubbard center
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom\s*%s.*\s*\(\s*l=%d.*\n"%(atmSp,ql),re.MULTILINE)          
                              

				eqOrbList = map(int, map(float,eqOrbRegex.findall(proj_lines)))#;print eqOrbList

				red_basis_for2e = [x - 1 for x in eqOrbList]
				#Get list of orbitals of type l for one atom of the species
				red_basis_2e = []
				red_basis_2e.append(red_basis_for2e[0])
				for i in range(1,len(red_basis_for2e)):
					if float(red_basis_for2e[i]) == float(red_basis_for2e[i-1])+1:
                                            red_basis_2e.append(red_basis_for2e[i])
					else:break

				#Create input file for respective species
				infnm = "_" + oneCalcID + "_acbn0_infile_%s.txt"%atmSp
				fout = file(os.path.join(subdir,infnm), 'w')
				S = "latvects = " + cellParaStr + "\n"
				fout.write(S)   
				S = "coords = " + atmPosStr + "\n"
				fout.write(S)
				S = "atlabels = " + atmLblsStr + "\n"
				fout.write(S)
				fout.write("nspin = %d\n" % nspin)
				fout.write("fpath = %s\n" % subdir)
				outfnm = "_" + oneCalcID + "_acbn0_outfile_%s.txt"%atmSp
				fout.write("outfile = %s\n"%outfnm)
				S = "reduced_basis_dm = " + str(red_basis).strip('[]') + "\n"
				fout.write(S)
				S = "reduced_basis_2e = " + str(red_basis_2e).strip('[]') + "\n"
				fout.write(S)
				fout.close()

				#Add filename to acbn0 run list
				inFileList.append(infnm)

			except Exception,e:
				AFLOWpi.run._fancy_error_log(e)

		return inFileList
			

	def run_acbn0(inputFiles):
			
		for infnm in inputFiles:
			
			cmd="python %s/acbn0.py %s > /dev/null"%(subdir,os.path.join(subdir,infnm))	
			print "Starting python acbn0.py %s\n"%(os.path.join(subdir,infnm))
			logging.info("Starting python acbn0.py %s\n"%(os.path.join(subdir,infnm)))
			try:
				os.system(cmd)
		
			except Exception,e:
				AFLOWpi.run._fancy_error_log(e)
			print "Finished python acbn0.py %s\n"%(os.path.join(subdir,infnm))
			logging.info("Finished python acbn0.py %s\n"%(os.path.join(subdir,infnm)))
	acbn0_inFileList = gen_input(oneCalcID,subdir,nspin)	
	run_acbn0(acbn0_inFileList)


def getU_frmACBN0out(oneCalc,ID,byAtom=False):

        oneCalcID = ID#oneCalc['_AFLOWPI_PREFIX_'][1:] 
        subdir = oneCalc['_AFLOWPI_FOLDER_'] 

	#Get species
	inputfile =oneCalc['_AFLOWPI_INPUT_']
        if byAtom==False:
            species=list(set(AFLOWpi.retr._getPosLabels(inputfile)))
            uvalLogName='_uValLog.log'
        else:
            splitInput = AFLOWpi.retr._splitInput(inputfile)
            species=[str(x) for x in range(1,int(splitInput['&system']['nat'])+1)]
            uvalLogName='_uValLog_byAtom.log'
        
	Uvals = {}
                
        for isp in species:
        	#Check for acbn0 output in the work directory
		try:
			acbn0_outFile = subdir + "/_" + oneCalcID + "_acbn0_outfile_" + isp + ".txt"
			if os.path.isfile(acbn0_outFile):
				#Get U value from acbn0 output
				try:
					lines = file(acbn0_outFile, 'r').read()
					acbn0_Uval = re.findall("U_eff\s*=\s*([-]*\d+.\d+)",lines)[0]

					Uvals[isp] = float(acbn0_Uval)
				except Exception,e:
					AFLOWpi.run._fancy_error_log(e)
					print e
			else:
				Uvals[isp] = 0.001

		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)

	try:
		if os.path.isfile(os.path.join(subdir,'%s%s' % (ID,uvalLogName))):
			with open(os.path.join(subdir,'%s%s' % (ID,uvalLogName)),'a') as uValLog:
				uValLog.write('%s\n' % Uvals)
		else:
			with open(os.path.join(subdir,'%s%s' % (ID,uvalLogName)),'w') as uValLog:
				uValLog.write('%s\n' % Uvals)
	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)

	return Uvals


def run(calcs,uThresh=0.001,nIters=20,mixing=0.70,kp_mult=1.5):

    temp_calcs=copy.deepcopy(calcs)
    for ID,oneCalc in calcs.iteritems():
	    uThreshDict={}
	    for key in oneCalc.keys():

                if re.search(r'_AFLOWPI_[A-Z][0-9]*_', key):
                    
                    uThreshDict[oneCalc[key]]=uThresh
            
            '''
            in the chance that the input file to the frame already has
            U vals for the species set the initial uvaldict to those 
            '''
            try:
                splitInput=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
                spec_order = splitInput['ATOMIC_SPECIES']['__content__'].split('\n')
                isolated_spec=[]
                for x in spec_order:
                    try:
                        isolated_spec.append(x.split()[0].strip())
                    except:
                        pass                  
                for i in range(len(isolated_spec)):
                    try:
                        initial_U = splitInput['&system']['Hubbard_U(%s)'%(i+1)]
                        uThreshDict[isolated_spec[i]] = float(initial_U)
                    except:
                        pass
            except Exception,e:
                print e

	    loopblock = '''

uValue = %s

oneCalc, newUvals = AFLOWpi.scfuj._run(__submitNodeName__,oneCalc,ID,config=configFile,mixing=%s,kp_mult=%s)
print "New U values ", str(newUvals).strip('{}')
logging.info('newUvals = {0}'.format(newUvals))
for key in uValue.keys():

    try:
        if abs(uValue[key]-newUvals[key]) > %s and oneCalc['_SCFUJ_LoopCount_'] != %s:
           logging.info("scfuj did not converge, starting next iteration")
           print "scfuj did not converge, starting next iteration"
           oneCalc['__status__']['Complete']=False
           AFLOWpi.prep._saveOneCalc(oneCalc,ID)
           AFLOWpi.run._submitJob(ID,oneCalc,__submitNodeName__)
           sys.exit(0)

        elif abs(uValue[key]-newUvals[key]) > %s and oneCalc['_SCFUJ_LoopCount_'] == %s:
           logging.info("Maximum no. of iterations reached. scfuj did not converge")
           print "Maximum no. of iterations reached. scfuj did not converge"
           oneCalc['__status__']['Complete']=False
           oneCalc['__status__']['Error']='SCFUJ'
           AFLOWpi.prep._saveOneCalc(oneCalc,ID)
           sys.exit(0)


    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
logging.info('Completed scfuj convergence for final U values = {0}'.format(newUvals))
print 'Completed scfuj convergence for final U values = ',str(newUvals).strip('{}')

AFLOWpi.scfuj._get_ham_xml(oneCalc,ID)
'''  % (uThreshDict,mixing,kp_mult,uThresh,nIters,uThresh,nIters-1)
            if AFLOWpi.prep._findInBlock(oneCalc,ID,block='RUN',string='uValue') == False:
                AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',loopblock)
            if AFLOWpi.prep._findInBlock(oneCalc,ID,block='IMPORT',string='import sys') == False:
                AFLOWpi.prep._addToBlock(oneCalc,ID,'IMPORT','import sys\n')
                        
    AFLOWpi.run._skeletonRun(calcs)




# def _assignLabels(oneCalc,ID,positions_only=False):

#     oneCalcOrig=copy.deepcopy(oneCalc)
#     oneCalc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],oneCalc['_AFLOWPI_PREFIX_'][1:])
    
#     inputFileString = oneCalc['_AFLOWPI_INPUT_']

#     availableLetters =  list(string.ascii_uppercase+string.digits)
#     availableFirstLetters =  list(string.ascii_uppercase)

#     countNames = [''.join(x) for x in list(itertools.product(availableLetters,availableLetters))]
#     inputCalc = AFLOWpi.retr._splitInput(inputFileString)
#     labels = AFLOWpi.retr._getPosLabels(inputFileString)
    
#     speciesString = inputCalc['ATOMIC_SPECIES']['__content__']
#     speciesStringSplit=[x for x in speciesString.split('\n') if len(x.strip())!=0]
#     speciesList = sorted(list(set(AFLOWpi.retr._getPosLabels(inputFileString))))
#     speciesMapping=OrderedDict()

#     for items in range(len(speciesList)):
#         speciesMapping[speciesList[items]]=[availableFirstLetters[items],0,[]]

    
#     labelMapping=OrderedDict()
#     positionString = inputCalc['ATOMIC_POSITIONS']['__content__']
#     positionStringSplit = [x for x in positionString.split('\n') if len(x.strip())!=0]

#     for position in range(len(positionStringSplit)):
#         posSplit=positionStringSplit[position].split()
#         species=posSplit[0].strip()
#         count = speciesMapping[species][1]
#         newLabel = speciesMapping[species][0]+countNames[count]
        
#         labelMapping[newLabel]=species
#         speciesMapping[species][2].append(newLabel)
#         positionStringSplit[position]=re.sub(species,newLabel,positionStringSplit[position])
#         speciesMapping[species][1]+=1



#     positionStringNew = '\n'.join(positionStringSplit)

#     inputCalc['ATOMIC_POSITIONS']['__content__']=positionStringNew

#     newSpeciesStringList=[]
#     for item1 in range(len(speciesStringSplit)):
#         splitSpecies = speciesStringSplit[item1].split()
#         species = splitSpecies[0].strip()
#         mapList = speciesMapping[species][2]

#         for items in mapList:
            
#             newSpeciesEntry = speciesStringSplit[item1].split()
#             newSpeciesEntry[0]=items
#             newSpeciesStringList.append(' '.join(newSpeciesEntry))

#     newSpeciesString = '\n'.join(newSpeciesStringList)
#     newNTYP = len(newSpeciesStringList)
#     inputCalc['ATOMIC_SPECIES']['__content__'] = newSpeciesString
#     inputCalc['&system']['ntyp'] = '%s' % newNTYP
#     newInputString=AFLOWpi.retr._joinInput(inputCalc)

#     oneCalc['_AFLOWPI_INPUT_'] = newInputString



#     oneCalc['__scfuj_label_mapping__']=labelMapping
#     oneCalcOrig['__scfuj_label_mapping__']=labelMapping

#     with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % oneCalc['_AFLOWPI_PREFIX_'][1:]),'w') as inputFile:
#         inputFile.write(oneCalc['_AFLOWPI_INPUT_'])
#     with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % ID),'w') as inputFile:
#         inputFile.write(oneCalc['_AFLOWPI_INPUT_'])

# #        AFLOWpi.prep._saveOneCalc(oneCalc,oneCalc['_AFLOWPI_PREFIX_'][1:])
#         AFLOWpi.prep._saveOneCalc(oneCalc,ID)

#     oneCalcOrig['_AFLOWPI_INPUT_']=newInputString
#     return oneCalcOrig


def _run(__submitNodeName__,oneCalc,ID,config=None,mixing=0.70,kp_mult=1.5):
	execPrefix = ''
	execPostfix = ''
        oneCalcID = ID


        def abortIFRuntimeError(subdir, ID):
            outfile = file(os.path.join(subdir, "%s.out"%ID)).read()
            errorList = re.findall(r'from (.*) : error #.*\n',outfile)
            if len(errorList) > 0:        
                logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
                print "Error in %s.out -- ABORTING ACBN0 LOOP"%ID                    
                raise SystemExit



        if '__runList__' not in oneCalc.keys():
            oneCalc['__runList__']=[]

            
	if config!=None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug('forced config %s' % config)
	else:
		try:
			config = AFLOWpi.prep._getConfigFile()
			AFLOWpi.prep._forceGlobalConfigFile(config)
		except Exception,e:
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

            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)


##################################################################################################################
            AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)


            oneCalc['__runList__'].append('scf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            

            nscf_calc,nscf_ID= nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_mult)	



        else:
            '''if we are restarting from a job killed from going walltime 
            try to load ID_nscf and if we can't then just make a new one'''
            try:
                nscf_ID='%s_nscf' % ID
                nscf_calc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],nscf_ID)                

                '''we have to make sure nscf step has the correct walltime and start time if it's a restart'''
                nscf_calc['__walltime_dict__']=oneCalc['__walltime_dict__']
            except Exception,e:
                try:
                    nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_mult)	

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)



##################################################################################################################
        if 'nscf' not in oneCalc['__runList__']:

            try:
                npool=AFLOWpi.retr._get_pool_num(nscf_calc,nscf_ID)        

                if npool!=1:
                    if len(re.findall(r'npool\s*(?:\d+)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool\s*(?:\d+)','npool %s'%npool,execPostfix)
                        logging.debug(execPostfixPrime)

            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)


            AFLOWpi.run._oneRun(__submitNodeName__,nscf_calc,nscf_ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)

            AFLOWpi.retr._writeEfermi(nscf_calc,nscf_ID)

            abortIFRuntimeError(subdir, nscf_ID)
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            oneCalc['__runList__'].append('nscf')
	
##################################################################################################################
        pdos_calc,pdos_ID = AFLOWpi.scfuj.projwfc(oneCalc,ID)


        if not re.match('northo',execPostfix) or not re.match('no',execPostfix):
            execPostfix+=' -northo 1'

        if 'pdos' not in oneCalc['__runList__']:
            pdosPath = os.path.join(AFLOWpi.prep._ConfigSectionMap('prep','engine_dir'),'projwfc.x')

            AFLOWpi.run._oneRun(__submitNodeName__,pdos_calc,pdos_ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='custom',executable='projwfc.x',execPath=pdosPath)

            AFLOWpi.prep._form_TB_dir(oneCalc,ID)
#############
            oneCalc['__runList__'].append('pdos')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            abortIFRuntimeError(subdir, pdos_ID)



        splitInput = AFLOWpi.retr._splitInput(nscf_calc['_AFLOWPI_INPUT_'])
        
#        nbnd = int(splitInput['&system']['nbnd'])


        
        if 'bands' not in oneCalc['__runList__']:
            want_dict = AFLOWpi.scfuj.WanT_bands(oneCalc,ID)




            for want_ID,want_calc in want_dict.iteritems():
#                AFLOWpi.run._oneRun(__submitNodeName__,want_calc,want_ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='custom',execPath='./want_bands.x' )
                AFLOWpi.run._oneRun(__submitNodeName__,want_calc,want_ID,execPrefix="",execPostfix='',engine='espresso',calcType='custom',execPath='./want_bands.x' )

            AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham_up.txt")
            AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham_down.txt")
            AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham.txt")
            AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp_up.txt")
            AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp_down.txt")
            AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp.txt")
            '''in case we're working in local scratch'''            
#           AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['ham','wan','space'])

            for want_ID,want_calc in want_dict.iteritems():
                abortIFRuntimeError(subdir, want_ID)


            oneCalc['__runList__'].append('bands')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)

            '''
            Will need to be filled in for the executable name with whatever wanT executable is called 
            and that executable needs to be moved to the calculation directory tree before this is called
            '''

############
##################################################################################################################

        #Get new U values from acbn0.py 
        try:
            acbn0(oneCalc, pdos_ID)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            raise SystemExit

#        acbn0(oneCalc, pdos_ID,byAtom=True)
#        getU_frmACBN0out(oneCalc,ID,byAtom=True)
        newUvals = getU_frmACBN0out(oneCalc,ID)

        #slight to help with oscillation
        try:
            old_U = oneCalc['_AFLOWPI_UVALS_']
        except:
            old_U = newUvals

#        mixing=0.70
        for spec,val in old_U.iteritems():
            newUvals[spec]=mixing*old_U[spec]+(1.0-mixing)*newUvals[spec]


        oneCalc['_AFLOWPI_UVALS_']=newUvals

        #Update Uvals
        oneCalc = updateUvals(oneCalc, newUvals,ID=ID)
	#update Uvals in _<ID>.py
	AFLOWpi.prep._modifyVarVal(oneCalc,ID,varName='uValue',value=newUvals)

        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

	a = ID + '.in'
        
	inputfile = oneCalc['_AFLOWPI_INPUT_']
        with  open(os.path.join(subdir,a),'w') as new_inputfile:
		new_inputfile.write(inputfile)

        '''make sure the input is set to 'from_scratch'''
        AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','restart_mode',"'from_scratch'")

	'''
	return the new uVals to the main script where the 
        conditional statement will decide to rerun with 
        the new U vals as input or keep going on with the 
        rest of the script.
	'''

        '''save one calc with the full list so next iteration it runs everything fresh'''
        oneCalc['__runList__']=[]


        try:
            oneCalc['_SCFUJ_LoopCount_'] += 1
        except:
            oneCalc['_SCFUJ_LoopCount_']  = 1

        oneCalc['__status__']['SCFUJ Iteration']=oneCalc['_SCFUJ_LoopCount_']

#        try:
#            '''check for uval oscillation and if it's happening end the loop'''
#            if checkOscillation(ID,new_oneCalc)==True:
#                new_oneCalc['__status__']['Error']='U Value Oscillation'
#                AFLOWpi.prep._saveOneCalc(new_OneCalc,ID)
#                raise SystemExit
#        except Exception,e:
#            AFLOWpi.run._fancy_error_log(e)
            


        '''remove all the extra .oneCalc and input files files generated to start fresh'''
        for stage in ['nscf','pdos','WanT_bands']:
            try:
                os.remove('./_%s_%s.oneCalc'%(ID,stage))
#                os.remove('./%s_%s.in'%(ID,stage))
            except:
                pass
            

        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        try:
            pass
#            oneCalc,ID = AFLOWpi.prep._oneUpdateStructs(oneCalc,ID,override_lock=True)

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


        try:
            AFLOWpi.prep._clean_want_bands(oneCalc,ID)
            AFLOWpi.plot._bands(oneCalc,ID,yLim=[-10,10],postfix='acbn0_TB',tight_banding=True)
        except:
            pass


        return oneCalc, Uvals


def _get_ham_xml(oneCalc,ID):
    tempdir=AFLOWpi.prep._get_tempdir()

    wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
    if os.path.isabs(wantdir) == False:
        configFileLocation = AFLOWpi.prep._getConfigFile()
        configFileLocation = os.path.dirname(configFileLocation)
        wantdir =  os.path.join(configFileLocation, wantdir)

    iotk_exe=os.path.join(wantdir,'iotk')
    command = iotk_exe
    command += ' --iotk-exe %s.x' % iotk_exe
    command += ' convert ./%s_TB_WanT.ham ./RHAM_%s_%s.xml' % (ID,AFLOWpi.retr._getStoicName(oneCalc),ID)
    os.system(command)

def checkOscillation(ID,oneCalc,uThresh=0.001):
    
    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_uValLog.log'% ID),'r') as uValLogFile:
        uValLogString=uValLogFile.read()
    iterationList = [x for x in uValLogString.split('\n') if len(x)!=0]
    iterationDictList=[]
    for fileString in iterationList:
            iterationDictList.append(ast.literal_eval(fileString))
    if len(iterationDictList) < 3:
        return False
    species = iterationDictList[0].keys()
    bySpecies=collections.OrderedDict()
    for i in range(len(iterationDictList)):
        for U in range(len(species)):
            uVal = iterationDictList[i][species[U]]
            try:
                bySpecies[species[U]].append(uVal)
            except:
                bySpecies[species[U]]=[uVal]
#    if bySpecies

    bySpeciesDiff=collections.OrderedDict()            
    for species,uVals in bySpecies.iteritems():
        try:
            bySpeciesDiff[species] = [j-i for i, j in zip(uVals[:-1], uVals[1:])]
        except:
            bySpeciesDiff[species]=[0]

    for species,diff in bySpeciesDiff.iteritems():
        if diff!=reversed(sorted(diff)):
            return True
    return False





############
##################################################################################################################

        
# def orderTest(diff):
#     for i in range(len(diff)-1):
#         if diff[i] < diff[i+1] and diff[i+1] > 0.1:
#             return False
#     return True

# def convTest(Diff=[],uThresh = 0.001):
#         for i in Diff:
#             if orderTest(i[-3:]) == False:
#                 return False, False
#             elif orderTest(i[-3:]) and i[-1] > uThresh:
#                 return True, False
#             elif orderTest(i[-3:])== False and i[-1] < uThresh:
#                 return False,True
#         return True, True
