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

def _get_PAO_orb_count(ID,oneCalc):
    tmpdir = AFLOWpi.prep._get_tempdir()
    dfxml  = os.path.join(tmpdir,"%s.save/data-file.xml"%oneCalc["_AFLOWPI_PREFIX_"])

    if  os.path.isfile(dfxml):
        with open(dfxml) as ofo:
            ofs = ofo.read()

        tot_orbs = int(re.findall("NUMBER_OF_ATOMIC.*\n\s*(\d+)\s*",ofs)[-1])
    else:
        dfxml  = os.path.join(tmpdir,"%s.save/data-file-schema.xml"%oneCalc["_AFLOWPI_PREFIX_"])
        with open(dfxml) as ofo:
            ofs = ofo.read()
        tot_orbs = int(re.findall(        '<num_of_atomic_wfc>(\d+)</num_of_atomic_wfc>',ofs)[-1])
    
    return tot_orbs

def _want_txt_to_bin(fpath,fname):

    fns=fname.split(".")
    bin_file = os.path.join(fpath,fns[0]+".npy")
    if not os.path.exists( os.path.join(fpath,fname)):
        return

    fin   = open(fpath+'/'+fname,"r")
    # fs = fin.read()
    # if len(re.findall('[*]',fs))!=0:
    #     fin.close()
    #     fs = re.sub('.*[*]+.*\n','0.0 0.0\n',fs)
    #     fin   = open(fpath+'/'+fname,"w")
    #     fin.write(fs)
    #     fin.close()
    #     fin   = open(fpath+'/'+fname,"r")
    try:
        ret=np.asarray(list(csv.reader(fin, delimiter=' ',skipinitialspace=True,
                                       quoting=csv.QUOTE_NONNUMERIC)),dtype=np.float32)
        fin.close()

        bin_data = ret[:,0]+1j*ret[:,1]

        bout   = open(bin_file,"wb")
        np.save(bout,bin_data)
        bout.close()
    except:
        return




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
        if PAOdir is None:
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
def scfprep(calcs,paodir=None,U_eff=True):
    output_calcs=copy.deepcopy(calcs)
    for ID,oneCalc in list(calcs.items()):
        AFLOWpi.prep._addToBlock(oneCalc,ID,'PREPROCESSING','oneCalc,ID = AFLOWpi.scfuj._oneScfprep(oneCalc,ID,U_eff=%s)'%U_eff) 
    return output_calcs

def _oneScfprep(oneCalc,ID,paodir=None,U_eff=True):
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

        if 'SCFUJ Iteration' in  list(oneCalc['__status__'].keys()):
            return oneCalc,ID

        configFileLocation = AFLOWpi.prep._getConfigFile()
        if not os.path.exists(configFileLocation):
            configFileLocation = oneCalc['_AFLOWPI_CONFIG_']
#            break
        
        output_oneCalc = copy.deepcopy(oneCalc)
        if paodir is None:
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
#        temp_species,_ = AFLOWpi.prep._resolve_AS_order(inputfile)
#        species = temp_species.keys()
#        species = re.findall("(.*).*UPF",inputfile)
#        species = [x.split()[0] for x in species]
        species=[x.split()[0] for x in AFLOWpi.retr._splitInput(inputfile)["ATOMIC_SPECIES"]["__content__"].split("\n") if len(x.strip()) !=0]
        splitInput = AFLOWpi.retr._splitInput(inputfile)


        Uvals = {}
        Jvals = {}
        uppered_system={}
        for k,v in list(splitInput['&system'].items()):
            uppered_system[k.upper()]=v


        for isp in range(len(species)):
            if 'HUBBARD_U(%s)' % (isp+1) in list(uppered_system.keys()):
                startingUval = float(uppered_system['HUBBARD_U(%s)' % (isp+1)])
                Uvals[species[isp]] = startingUval
            else:
                Uvals[species[isp]] = 0.001

            if U_eff==False:
                if 'HUBBARD_J0(%s)' % (isp+1) in list(uppered_system.keys()):
                    startingJval = float(uppered_system['HUBBARD_J0(%s)' % (isp+1)])
                    Jvals[species[isp]] = startingJval
                else:
                    Jvals[species[isp]] = 0.001 



        AFLOWpi.prep._modifyVarVal(oneCalc,ID,varName='uValue',value=Uvals)
        AFLOWpi.prep._modifyVarVal(oneCalc,ID,varName='jValue',value=Jvals)
        new_calc = updateUvals(oneCalc,Uvals,Jvals,ID=ID,U_eff=U_eff)

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

def updateUvals(oneCalc, Uvals,Jvals,ID=None,U_eff=True):
        """
        Modify scf input file to do a lda+u calculation.
        
        Arguments:
         - oneCalc      -- Dictionary of one calculation 
         - Uvals        -- Dictionary of Uvals

        """
        try:
                inputfile = oneCalc['_AFLOWPI_INPUT_']
                #Get species
#                species=list(set(AFLOWpi.retr._getPosLabels(inputfile)))
#               species = re.findall("(\w+).*UPF",inputfile)

                species=[x.split()[0] for x in AFLOWpi.retr._splitInput(inputfile)["ATOMIC_SPECIES"]["__content__"].split("\n") if len(x.strip()) !=0]
#                species = re.findall("(.*).*UPF",inputfile)
#                species = [x.split()[0] for x in species]

#                AFLOWpi.prep._get_order_from_atomic_species(inputfile)
                inputDict = AFLOWpi.retr._splitInput(inputfile)
                inputDict['&system']['lda_plus_u']='.TRUE.'

                if "noncolin" in list(inputDict['&system'].keys()):
                    if inputDict['&system']["noncolin"]==".true.":
                        inputDict['&system']["lda_plus_u_kind"]=1
                    else:
                        inputDict['&system']["lda_plus_u_kind"]=0
                else:
                    inputDict['&system']["lda_plus_u_kind"]=0

                inputDict['&system']['lda_plus_u']

                for isp in range(len(species)):
                    hub_entry = 'Hubbard_U(%s)'% str(isp+1)
                    inputDict['&system'][hub_entry] = Uvals[species[isp]]
                    if U_eff==False:
                        hubj_entry = 'Hubbard_J0(%s)'% str(isp+1)
                        inputDict['&system'][hubj_entry] = Jvals[species[isp]]

                #Update inputfile
                oneCalc['_AFLOWPI_INPUT_']=AFLOWpi.retr._joinInput(inputDict)

        except Exception as e:
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
         - paodir       - path of pseudoatomic orbital basis set
        """


        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % ID),'w') as scfujInfileObj:
            scfujInfileObj.write(oneCalc['_AFLOWPI_INPUT_'])

        '''save the calc just in case so it's updated witht he scfuj U value info in the input file as well as _AFLOWPI_INPUT_'''



        espressodir = AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if os.path.isabs(espressodir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            espressodir =  os.path.join(configFileLocation, espressodir)

        pdosExec = 'projwfc.x'
        pdosExec = os.path.join(espressodir,pdosExec)
        if os.path.exists(pdosExec)==False:
            print('ProjectWFC executable not found. Check your config file to make sure engine_dir path is correct and that projwfc.x is in that directory..Exiting')
            logging.error('ProjectWFC executable not found. Check your config file to make sure engine_dir path is correct and that projwfc.x is in that directory..Exiting')
            raise SystemExit
        scfExec = 'pw.x'
        scfExec = os.path.join(espressodir,scfExec)
        if os.path.exists(scfExec)==False:
            print('SCF executable not found. Check your config file to make sure engine_dir path is correct and that pw.x is in that directory..Exiting')
            logging.error('SCF executable not found. Check your config file to make sure engine_dir path is correct and that pw.x is in that directory..Exiting')
            raise SystemExit
            




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
            if paodir is None:
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
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
            print('Can not find PAO file. Exiting')
            logging.error('Can not find PAO file. Exiting')
            raise SystemExit




        #Copy PAO files 

        try:
            work_dir = oneCalc['_AFLOWPI_FOLDER_'] 
            for key in list(oneCalc.keys()):
                    v = oneCalc[key]
                    if re.search(r'_AFLOWPI_[A-Z][0-9]*_', key):
                            vp = AFLOWpi.scfuj._getPAOfilename(v.strip('0123456789'),paodir)            


                            try:
                                    a = os.path.join(paodir,vp)
                            except AttributeError:
                                if v=='!':
                                    continue
                                logging.debug('Cannot find correct PAO files in %s ...Exiting' % paodir)
                                print(('cannot find correct PAO files in %s ...Exiting' % paodir))
                                raise SystemExit

                            newPAOFileNm = os.path.join(work_dir,v.strip('0123456789')+"_basis.py")
                            print(('Copying '+a+' to '+ newPAOFileNm+'\n'))
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

        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


        '''save the inital list of things to exec so first iteration it will run through'''


        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        return oneCalc

def nscf_nosym_noinv(oneCalc,ID=None,kpFactor=1.50,unoccupied_states=False,band_factor=1.25,tetra_nscf=False,wsyminv=False):
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

                        '''do nbnd==number of projectors'''
                        nbnd = AFLOWpi.scfuj._get_PAO_orb_count(ID,oneCalc)
                        splitInput['&system']['nbnd']=int(float(nbnd)*band_factor)

                        '''checks to see if nbnd has already been added to the file - in any case add/replace nbnd'''


                        '''Add nosym = .true. and noinv = .true. to file '''

                        try:
                            if "occupations" in list(splitInput['&system'].keys()) and tetra_nscf:
                                splitInput['&system']['occupations']='"tetrahedra"'
                                try:
                                    del splitInput['&system']['smearing']
                                except: pass
                                try:
                                    del splitInput['&system']['degauss']
                                except: pass

                            splitInput['&electrons']['conv_thr']='1.0D-6'

#                            if not wsyminv:
#                                splitInput['&system']['nosym']='.True.'
#                                splitInput['&system']['noinv']='.True.'

                            splitInput['&control']['verbosity']='"high"'
#                            splitInput['&control']['wf_collect']='.TRUE.'
                            splitInput['&control']['calculation']='"nscf"'
                            inputfile=AFLOWpi.retr._joinInput(splitInput)
                            '''writes an input for band_plot.x to process the correct number of bands calculated'''
                        
                            
                        except Exception as e:
                                AFLOWpi.run._fancy_error_log(e)
                                

                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)
                        print('SCF outputfile not found: nbnd is default')
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

#                        before_kpf=reduce(lambda x, y: x*y, scfKPointSplit )
#                        scaling_kpf=before_kpf/150.0
#                        if scaling_kpf>1.0:
#                            kpFactor=1
                        
                            
#                        kpFactor=2.0    

                        for kpoint in range(len(scfKPointSplit)-3):
                            try:
                                kpFactor[0]
                                scfKPointSplit[kpoint] = str(int(math.ceil(scfKPointSplit[kpoint]\
                                                                           *kpFactor[kpoint])))
                            except:
                                scfKPointSplit[kpoint] = str(int(math.ceil(scfKPointSplit[kpoint]*kpFactor)))
                        for kpoint in range(3,len(scfKPointSplit)):
                                scfKPointSplit[kpoint] = '0'

                        newKPointString = ' '.join(scfKPointSplit)

                        splitInput['K_POINTS']['__content__']=newKPointString
                        splitInput['K_POINTS']['__modifier__']='{automatic}'
                        inputfile = AFLOWpi.retr._joinInput(splitInput)

                except Exception as e:
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


def projwfc(oneCalc,ID=None,paw=False,ovp=False):
        '''
        Run projwfc on each calculation

        Arguments:
         - oneCalc -- dictionary of a single calculation

        '''

        temp_dir= AFLOWpi.prep._get_tempdir()
        calc_copy=copy.deepcopy(oneCalc)
        output_calc = {}
        nscf_ID=ID+'_nscf'
        Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
        eShift=float(Efermi)+4.0
        if ovp:
            ovp_str=".TRUE."
        else:
            ovp_str=".FALSE."

        prefix = oneCalc['_AFLOWPI_PREFIX_']
        try:
                subdir = oneCalc['_AFLOWPI_FOLDER_']
                try:
                    prefix = AFLOWpi.retr._prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])
                except Exception as e:
                    prefix = oneCalc['_AFLOWPI_PREFIX_']
                inputfile = """&PROJWFC
  prefix='%s'
  filpdos='./%s'
  outdir='%s'
  lwrite_overlaps=%s
  lbinary_data  = .FALSE.
/
"""%(prefix,ID,temp_dir,ovp_str)

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
        si=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
        try:
            nspin = int(si['&system']['nspin'])
            return nspin
        except Exception as e:
        #    if "noncolin" in si['&system'].keys():
        #        if "magnetization" in " ".join(si['&system'].keys()):
        #            return 2
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

            except Exception as e:
                pass


        

def evCurveMinimize(calcs,config=None,pThresh=10.0,final_minimization = 'vc-relax'):

        for ID,oneCalc in list(calcs.items()):
                execString = '''
newOneCalc = AFLOWpi.scfuj. _oneMinimizeCalcs(oneCalc, ID, pThresh=%f)
AFLOWpi.prep._saveOneCalc(oneCalc,ID)
        '''%(pThresh)
                AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN', execString)
        AFLOWpi.run._skeletonRun(calcs)
        #if final_minimization is not None:
#               calcs = AFLOWpi.prep.changeCalcs(calcs, 'calculation', final_minimization)
#               AFLOWpi.run.scf(calcs)
#               calcs = AFLOWpi.prep.changeCalcs(calcs, 'calculation', 'scf')
#       AFLOWpi.run.scf(calcs)

        return calcs

def _oneMinimizeCalcs(oneCalc,ID,config=None,pThresh=10.0):
        '''
                Get equilibrium volume using evfit to Murnaghan's EOS with 5 volumes in +/-20% of input volume
        '''



        __submitNodeName__ =  __main__. __submitNodeName__
        execPrefix = ''
        execPostfix = ''
        if config is not None:
                AFLOWpi.prep._forceGlobalConfigFile(config)
                logging.debug('forced config %s' % config)
        else:
                try:
                        config = AFLOWpi.prep._getConfigFile()
                        AFLOWpi.prep._forceGlobalConfigFile(config)
                except Exception as e:
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
                                with open(outFile) as ifo:
                                    outFileString = ifo.read()
                                energyRegex = re.compile(r'(?:(?:(?:(?:\!\s+)total)|(?:Final)) en\w+\s*=\s+(.+?)Ry)',re.MULTILINE)
                                finalEnergy=energyRegex.findall(outFileString)
                                stressRegex = re.compile(r'P=\s*(.*)\n')
                                finalStress = stressRegex.findall(outFileString)
                                

                                if len(finalEnergy):
                                        energy  = float(finalEnergy[-1]); engList.append(energy)
                                        stress = float(finalStress[0]); stressList.append(stress)
                                        logging.info("E-V curve point: %f a.u. %f Ry" %(i, energy))
                                        print(("E-V curve point: %f a.u. %f Ry" %(i, energy)))  
                                else:
                                        logging.error("ERROR: E-V curve at celldm(1) = %f FAILED"%i)
                                        raise SystemExit

                        except Exception as e:
                                AFLOWpi.run._fancy_error_log(e)

                return list(zip(stressList, engList, alatList))

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
                        with open(evinFile,'a') as ifo:
                            fout = ifo.read()
                        dataList.sort(key=lambda tup: tup[2])
                        #Write data to file
                        for i in dataList:
                                fout.write("%f\t %f\n"%(i[2], i[1]))
                        fout.close()

                        engineDir=AFLOWpi.prep._ConfigSectionMap("prep","engine_dir")
                        evfitPath=os.path.join(engineDir,'ev.x')
                        evfitString = "cat<<! | %s \nau\nsc\n4\n%s\n%s\n!\n"%(evfitPath,evinFile,evoutFile)
                        logging.info("Starting fitting E-V curve with Murnaghan's EOS")
                        print("Starting Fitting E-V curve with Murnaghan's EOS")
                        os.system(evfitString)
                        with open(evoutFile,'r') as ifo:
                            evoutFileString = ifo.read()
                        logging.info("Finished fitting E-V curve with Murnaghan's EOS")
                        print("Finished Fitting E-V curve with Murnaghan's EOS")
                        
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


                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)

                return newCalc
        else:
                print(("Minimization for ibrav = %d not implemented" %ibrav))
                logging.error("Minimization for ibrav = %d not implemented"%ibrav)
                raise SystemExit
def acbn0(oneCalc,projCalcID,byAtom=False):
        '''
        



        '''

        oneCalcID = '_'.join(projCalcID.split('_')[:-1])#oneCalc['_AFLOWPI_PREFIX_'][1:]        

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
                        a,cellParaMatrix=AFLOWpi.retr._getCellParams(oneCalc,oneCalcID)
                        ALAT = float(AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])["&system"]["celldm(1)"])


                        """THINK OF A BETTER WAY TO CHECK THIS"""



                        l=(cellParaMatrix*a).tolist()
                        cellParaStr = ""
                        """THINK OF A BETTER WAY TO CHECK THIS"""
                        for i in range(3):
                                cellParaStr += str(l[i]).strip('[]')
                                if i != 2:cellParaStr += ' ,'

                        #Get atomic positions in cartesian coordinates, in a single string in pattern x1,y1,z1, x2, y2, z2, ..., xn, yn, zn
                        scfOutput = '%s.out' % oneCalcID
                        with open(os.path.join(subdir,scfOutput),'r') as ifo:
                            lines = ifo.read()

                        splitInput = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
                        atmPosRegex = re.compile(r"positions \(alat units\)\n((?:.*\w*\s*tau\(.*\)\s=.*\(.*\)\n)+)",re.MULTILINE) 
                        FROM_INPUT=False
                        try:
                            positions = AFLOWpi.retr.getPositionsFromOutput(oneCalc,oneCalcID)
                            positions = AFLOWpi.retr._cellMatrixToString(positions)
                            
                            if positions=="":
                                FROM_INPUT=True
                                splitInput = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
                                positions= splitInput["ATOMIC_POSITIONS"]["__content__"]
                        except:
                            FROM_INPUT=True
                            splitInput = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
                            positions= splitInput["ATOMIC_POSITIONS"]["__content__"]

                        if len(positions.strip())==0:
                            positions= splitInput["ATOMIC_POSITIONS"]["__content__"]


                        lines1 = atmPosRegex.findall(lines)[0]
                        atmPosRegex1 = re.compile(r".*=.*\((.*)\)\n+",re.MULTILINE)

                        atmPos = atmPosRegex1.findall(lines1)


                        try:
                            positions=np.asarray([[float(i.split()[1]),float(i.split()[2]),float(i.split()[3])] for i in positions.split("\n") if len(i.strip())!=0])
                        except:
                            positions=np.asarray([[float(i.split()[0]),float(i.split()[1]),float(i.split()[2])] for i in positions.split("\n") if len(i.strip())!=0])



                        positions=AFLOWpi.retr._convertCartesian(positions,cellParaMatrix,scaleFactor=1.0)
                        #in bohr


                        positions*=a
                       

                        atmPos=AFLOWpi.retr._cellMatrixToString(positions).split("\n")
                        
                        atmPosList = []

                        for i in atmPos:atmPosList.append(list(map(float,i.split())))
                        atmPosStr = ""

                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)
                        raise SystemExit


                for i in range(len(atmPosList)):
                        try:
                #Convert fractional atomic coordinates to cartesian coordinates using the lattice vectors (cell parameters)
#                               atmPosStr += str(list(np.array(atmPosList[i])*a)).strip('[]')
                                atmPosStr += str(list(np.array(atmPosList[i]))).strip('[]')
                                if i != len(atmPosList)-1:
                                        atmPosStr += ' ,'
                        except Exception as e:
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


                                #Get list of orbitals
                                projOut = projCalcID + '.out'
                                with open(os.path.join(subdir,projOut), 'r') as ifo:
                                    proj_lines = ifo.read()

                                inFileList = []

                        except Exception as e:
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
#                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*\(%s\d*\s*\).*\(l=%d.*\)\n"%(atmSp.strip('0123456789'),ql),re.MULTILINE)
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*\(%s\d*\s*\).*\(l=%d.*\)\n"%(atmSp,ql),re.MULTILINE)


                                    
                                    eqOrbList = list(map(int, list(map(float, eqOrbRegex.findall(proj_lines)))))
                                    red_basis = [x - 1 for x in eqOrbList]

                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*\(%s\s*\).*\(l=%d.*\)\n"%(atmSp,ql),re.MULTILINE)
                                #Get ones relevant for hubbard center


                                else:

                                    getSpecByAtomRegex = re.compile(r"state #\s*(?:\d*): atom\s*%s\s*\(([a-zA-Z]+)"%atmSp)
                                    speciesFromNum = getSpecByAtomRegex.findall(proj_lines)[-1]
                                    ql = get_orbital(speciesFromNum.strip('0123456789'))

                                    #Get list of all orbitals of type ql of the same atom
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*(%s\s*).*l=%d.*\n"%(speciesFromNum,ql),re.MULTILINE)  
#                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom.*%s.*l=%d.*\n"%(speciesFromNum.strip('0123456789'),ql),re.MULTILINE)          


                                    eqOrbList = list(map(int, list(map(float,eqOrbRegex.findall(proj_lines)))))
                                    red_basis = [x - 1 for x in eqOrbList]

                                #Get ones relevant for hubbard center
                                    eqOrbRegex = re.compile(r"state #\s*(\d*): atom\s*%s.*\s*\(\s*l=%d.*\n"%(atmSp,ql),re.MULTILINE)          
                              

                                eqOrbList = list(map(int, list(map(float,eqOrbRegex.findall(proj_lines)))))#;print eqOrbList
                            
                                red_basis_for2e = [x - 1 for x in eqOrbList]
                                #Get list of orbitals of type l for one atom of the species
                                red_basis_2e = []
                                red_basis_2e.append(red_basis_for2e[0])
                                for i in range(1,len(red_basis_for2e)):
                                        if float(red_basis_for2e[i]) == float(red_basis_for2e[i-1])+1:
                                            red_basis_2e.append(red_basis_for2e[i])
                                        else:break

                                #Create input file for respective species
                                infnm = oneCalcID + "_acbn0_%s.in"%atmSp
                                fout = open(os.path.join(subdir,infnm), 'w')
                                S = "latvects = " + cellParaStr + "\n"
                                fout.write(S)   
                                S = "coords = " + atmPosStr + "\n"
                                fout.write(S)
                                S = "atlabels = " + atmLblsStr + "\n"
                                fout.write(S)
                                fout.write("nspin = %d\n" % nspin)
                                fout.write("fpath = %s\n" % subdir)
                                outfnm = oneCalcID + "_acbn0_%s.out"%atmSp
                                fout.write("outfile = %s\n"%outfnm)
                                S = "reduced_basis_dm = " + str(red_basis).strip('[]') + "\n"
                                fout.write(S)
                                S = "reduced_basis_2e = " + str(red_basis_2e).strip('[]') + "\n"
                                fout.write(S)
                                fout.close()

                                #Add filename to acbn0 run list
                                inFileList.append(infnm)

                        except Exception as e:
                                AFLOWpi.run._fancy_error_log(e)

                return inFileList
                        

        def run_acbn0(inputFiles):
                        
                for infnm in inputFiles:
                        
                        cmd="python %s/acbn0.py %s > /dev/null"%(subdir,os.path.join(subdir,infnm))     
                        print(("Starting python acbn0.py %s\n"%(os.path.join(subdir,infnm))))
                        logging.info("Starting python acbn0.py %s\n"%(os.path.join(subdir,infnm)))
                        try:
                                os.system(cmd)
                
                        except Exception as e:
                                AFLOWpi.run._fancy_error_log(e)
                        print(("Finished python acbn0.py %s\n"%(os.path.join(subdir,infnm))))
                        logging.info("Finished python acbn0.py %s\n"%(os.path.join(subdir,infnm)))
        acbn0_inFileList = gen_input(oneCalcID,subdir,nspin)    
        run_acbn0(acbn0_inFileList)


def getU_frmACBN0out(oneCalc,ID,byAtom=False,U_eff=True):

        oneCalcID = ID#oneCalc['_AFLOWPI_PREFIX_'][1:] 
        subdir = oneCalc['_AFLOWPI_FOLDER_'] 

        #Get species
        inputfile =oneCalc['_AFLOWPI_INPUT_']
        if byAtom==False:
#            species=list(set(AFLOWpi.retr._getPosLabels(inputfile)))
            temp_species,_ = AFLOWpi.prep._resolve_AS_order(inputfile)
            species = list(temp_species.keys())
            uvalLogName='_uValLog.log'
        else:
            splitInput = AFLOWpi.retr._splitInput(inputfile)
            species=[str(x) for x in range(1,int(splitInput['&system']['nat'])+1)]
            uvalLogName='_uValLog_byAtom.log'
        
        Uvals = {}
        Jvals = {}
                                    


        for isp in species:
                #Check for acbn0 output in the work directory
                try:
                        acbn0_outFile = subdir + "/" + oneCalcID + "_acbn0_" + isp + ".out"
                        if os.path.isfile(acbn0_outFile):
                                #Get U value from acbn0 output
                                try:
                                        with open(acbn0_outFile, 'r') as ifo:
                                            lines = ifo.read()
                                        acbn0_Uval = re.findall("U_eff\s*=\s*([-]*\d+.\d+)",lines)[0]
                                        acbn0_Jval = 0.0
                                        if not U_eff:
                                            acbn0_Uval = re.findall("Parameter U=\s*([-]*\d+.\d+)",lines)[0]
                                            try:
                                                acbn0_Jval = re.findall("Parameter J=\s*([-]*\d+.\d+)",lines)[0]
                                            except:
                                                acbn0_Jval="0.0"
                                        Uvals[isp] = float(acbn0_Uval)
                                        Jvals[isp] = float(acbn0_Jval)
                                except Exception as e:
                                        AFLOWpi.run._fancy_error_log(e)

                        else:
                                Uvals[isp] = 0.001
                                Jvals[isp] = 0.001

                except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)

        try:
                if os.path.isfile(os.path.join(subdir,'%s%s' % (ID,uvalLogName))):
                        with open(os.path.join(subdir,'%s%s' % (ID,uvalLogName)),'a') as uValLog:
                                uValLog.write('%s\n' % Uvals)
                else:
                        with open(os.path.join(subdir,'%s%s' % (ID,uvalLogName)),'w') as uValLog:
                                uValLog.write('%s\n' % Uvals)
        except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

        return Uvals,Jvals


def run(calcs,uThresh=0.001,nIters=20,mixing=0.10,kp_mult=1.6,U_eff=False):

    temp_calcs=copy.deepcopy(calcs)
    for ID,oneCalc in list(calcs.items()):
            uThreshDict={}
            jThreshDict={}
            for key in list(oneCalc.keys()):

                if re.search(r'_AFLOWPI_[A-Z][0-9]*_', key):
                    
                    uThreshDict[oneCalc[key]]=uThresh
                    if U_eff==False:
                        jThreshDict[oneCalc[key]]=uThresh/2.0
                    
            
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
                        if U_eff==False:
                            initial_J = splitInput['&system']['Hubbard_J0(%s)'%(i+1)]
                            jThreshDict[isolated_spec[i]] = float(initial_J)
                    except:
                        pass
            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

            loopblock = '''

uValue = %s

oneCalc, newUvals = AFLOWpi.scfuj._run(__submitNodeName__,oneCalc,ID,config=configFile,mixing=%s,kp_mult=%s,U_eff=%s)
print("New U values ", str(newUvals).strip('{}'))
logging.info('newUvals = {0}'.format(newUvals))
for key in uValue.keys():

    try:
        if abs(uValue[key]-newUvals[key]) > %s and oneCalc['_SCFUJ_LoopCount_'] != %s:
           logging.info("scfuj did not converge, starting next iteration")
           print("scfuj did not converge, starting next iteration")
           oneCalc['__status__']['Complete']=False
           AFLOWpi.prep._saveOneCalc(oneCalc,ID)
           AFLOWpi.run._submitJob(ID,oneCalc,__submitNodeName__)
           sys.exit(0)

        elif abs(uValue[key]-newUvals[key]) > %s and oneCalc['_SCFUJ_LoopCount_'] == %s:
           logging.info("Maximum no. of iterations reached. scfuj did not converge")
           print("Maximum no. of iterations reached. scfuj did not converge")
           oneCalc['__status__']['Complete']=False
           oneCalc['__status__']['Error']='SCFUJ'
           AFLOWpi.prep._saveOneCalc(oneCalc,ID)
           sys.exit(0)


    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
logging.info('Completed scfuj convergence for final U values = {0}'.format(newUvals))
print('Completed scfuj convergence for final U values = ',str(newUvals).strip('{}'))


'''  % (uThreshDict,mixing,kp_mult,U_eff,uThresh,nIters,uThresh,nIters-1)
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


def _run(__submitNodeName__,oneCalc,ID,config=None,mixing=0.10,kp_mult=1.6,U_eff=True):
        execPrefix = ''
        execPostfix = ''
        oneCalcID = ID


        def abortIFRuntimeError(subdir, ID):
            outfile = file(os.path.join(subdir, "%s.out"%ID)).read()
            errorList = re.findall(r'from (.*) : error #.*\n',outfile)
            if len(errorList) > 0:        
                logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
                print(("Error in %s.out -- ABORTING ACBN0 LOOP"%ID))                    
                raise SystemExit



        if '__runList__' not in list(oneCalc.keys()):
            oneCalc['__runList__']=[]

            
        if config is not None:
                AFLOWpi.prep._forceGlobalConfigFile(config)
                logging.debug('forced config %s' % config)
        else:
                try:
                        config = AFLOWpi.prep._getConfigFile()
                        AFLOWpi.prep._forceGlobalConfigFile(config)
                except Exception as e:
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

        nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,kpFactor=kp_mult,band_factor=1.0,wsyminv=False)    
##################################################################################################################
        
##################################################################################################################
        pdos_calc,pdos_ID = AFLOWpi.scfuj.projwfc(oneCalc,ID,ovp=True)


        if not re.match('northo',execPostfix) or not re.match('no',execPostfix):
            execPostfix+=' -northo 1'


        splitInput = AFLOWpi.retr._splitInput(nscf_calc['_AFLOWPI_INPUT_'])
        AFLOWpi.prep._run_tb_ham_prep(__submitNodeName__,oneCalc,ID,kp_factor=kp_mult,cond=0,ovp=True,band_factor=1.0,wsyminv=True)

        AFLOWpi.prep._from_local_scratch(oneCalc,ID,ext_list=['.save'])
        AFLOWpi.scfuj._add_paopy_header(oneCalc,ID,shift_type=1,shift=1.0,thresh=0.90,tb_kp_mult=1.0,acbn0=True,ovp=True,smearing='gauss')
        AFLOWpi.scfuj._run_paopy(oneCalc,ID,acbn0=True)

        AFLOWpi.prep._saveOneCalc(oneCalc,ID)


        '''
            Will need to be filled in for the executable name with whatever wanT executable is called 
            and that executable needs to be moved to the calculation directory tree before this is called
        '''
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham_up.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham_down.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kham.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp_up.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp_down.txt")
        AFLOWpi.scfuj._want_txt_to_bin(oneCalc['_AFLOWPI_FOLDER_'],"kovp.txt")

############
##################################################################################################################

        #Get new U values from acbn0.py 
        try:
            acbn0(oneCalc, pdos_ID)
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
            raise SystemExit



        newUvals,newJvals = getU_frmACBN0out(oneCalc,ID,U_eff=U_eff)
        Uvals=newUvals
        Jvals=newJvals

        #slight to help with oscillation
        try:
            old_U = oneCalc['_AFLOWPI_UVALS_']
        except:
            old_U = newUvals


        for spec,val in list(old_U.items()):
            newUvals[spec]=mixing*old_U[spec]+(1.0-mixing)*newUvals[spec]


        oneCalc['_AFLOWPI_UVALS_']=newUvals
        oneCalc['_AFLOWPI_JVALS_']=newJvals

        #Update Uvals
        oneCalc = updateUvals(oneCalc, newUvals,newJvals,ID=ID,U_eff=U_eff)
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
#        except Exception as e:
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

        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


        try:
#            AFLOWpi.prep._clean_want_bands(oneCalc,ID)
            AFLOWpi.plot._bands(oneCalc,ID,yLim=[-15,5],postfix='acbn0_TB',tight_banding=True)
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
    try:
        os.system(command)
    except:
        pass
    command = iotk_exe
    command += ' --iotk-exe %s.x' % iotk_exe
    command += ' convert ./%s_TB_WanT_up.ham ./RHAM_UP_%s_%s.xml' % (ID,AFLOWpi.retr._getStoicName(oneCalc),ID)
    try:
        os.system(command)
    except:
        pass
    command = iotk_exe
    command += ' --iotk-exe %s.x' % iotk_exe
    command += ' convert ./%s_TB_WanT_dn.ham ./RHAM_DOWN_%s_%s.xml' % (ID,AFLOWpi.retr._getStoicName(oneCalc),ID)
    try:
        os.system(command)
    except:
        pass


def checkOscillation(ID,oneCalc,uThresh=0.001):
    
    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_uValLog.log'% ID),'r') as uValLogFile:
        uValLogString=uValLogFile.read()
    iterationList = [x for x in uValLogString.split('\n') if len(x)!=0]
    iterationDictList=[]
    for fileString in iterationList:
            iterationDictList.append(ast.literal_eval(fileString))
    if len(iterationDictList) < 3:
        return False
    species = list(iterationDictList[0].keys())
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
    for species,uVals in list(bySpecies.items()):
        try:
            bySpeciesDiff[species] = [j-i for i, j in zip(uVals[:-1], uVals[1:])]
        except:
            bySpeciesDiff[species]=[0]

    for species,diff in list(bySpeciesDiff.items()):
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
