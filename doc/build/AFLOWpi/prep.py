import AFLOWpi

def __get_step(oneCalc,ID,step_type=None,last=True):

    if step_type==None:
        return 0
    try:
        workflow = oneCalc['_AFLOWPI_WORKFLOW_']
    except:
        return 0


    try:

        chain_index=oneCalc['__chain_index__']


        occ_list=[]
        current_step_type=workflow[chain_index-1]
        for i in reversed(range(0,chain_index)):

            if workflow[i]==step_type:


                occ_list.append(i+1)
                    
        if len(occ_list)!=0:
            return occ_list
        else:
            return 0


    except Exception,e:
        print e
        AFLOWpi.run._fancy_error_log(e)
        return 0
        


import numpy





def __return_ID(oneCalc,ID,step_type=None,last=True,straight=False):
    index =__get_step(oneCalc,ID,step_type=step_type,last=last)
    prefix = oneCalc['_AFLOWPI_PREFIX_'][1:]

    prefix_first = prefix.split('_')[0]



    if type(index)==type([1,2,3]):



        index = numpy.asarray(index)

        splits=numpy.split(index, numpy.where(numpy.diff(index) != -1)[0]+1)
        if last==True:
            chain_ind_list=splits[0].tolist()
            if straight==True:
                return ['%s_%02d'%(prefix_first,i) for i in chain_ind_list]
            else:
                return ['%s_%02d'%(prefix_first,i) for i in chain_ind_list][0]

        else:
            if straight==True:
                return ['%s_%02d'%(prefix_first,i) for i in index]
            else:
                return ['%s_%02d'%(prefix_first,i) for i in index[0]]

import AFLOWpi
import os
import __main__
import numpy
import time
import ConfigParser
import glob 
import sys 
import re

def __one_phon(oneCalc,ID,subset_name='SUBSET'):
    return glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/FD_PHONON'+'/displaced*.in')

def __one_test_build(oneCalc,ID,build_command,subset_name='SUBSET',merge_oneCalc=True,keep_name=False,config=None,clean_input=True):
    if config==None:
        config=oneCalc['_AFLOWPI_CONFIG_']

    intoInit={'PROJECT':subset_name,'SET':'','workdir':oneCalc['_AFLOWPI_FOLDER_'],'config':config}
    fake_session_keys = AFLOWpi.prep.init(**intoInit)

    exec('input_strings=%s'%build_command)
    if merge_oneCalc==True:
        varied_calcs = AFLOWpi.prep.calcFromFile(fake_session_keys,input_strings,reffile=oneCalc['_AFLOWPI_INPUT_'],workdir=oneCalc['_AFLOWPI_FOLDER_'],keep_name=keep_name,clean_input=clean_input)
    else:
        varied_calcs = AFLOWpi.prep.calcFromFile(fake_session_keys,input_strings,workdir=oneCalc['_AFLOWPI_FOLDER_'],keep_name=keep_name,clean_input=clean_input)

    return varied_calcs
###############################################################################################################

###############################################################################################################

def prep_split_step(calcs,subset_creator,subset_tasks=[],mult_jobs=False,substep_name='SUBSET',keep_file_names=False,clean_input=True):



#####################################################################
        AFLOWpi.run.__skeletonRun(calcs) 
	for ID,oneCalc in calcs.iteritems():

		oneCalc['__splitCounter__']=0

		execString='''if oneCalc['__execCounter__']<=%s:
    ''' % oneCalc['__execCounterBkgrd__']
		execString+='''
     oneCalc,ID = AFLOWpi.prep.construct_and_run(__submitNodeName__,oneCalc,ID,build_command=%s,subset_tasks=%s,mult_jobs=%s,subset_name='%s',keep_file_names=%s,clean_input=%s)

''' % (repr(subset_creator),repr(subset_tasks),mult_jobs,substep_name,keep_file_names,clean_input)  
#''' % (subset_creator,exit_command,repr(subset_tasks),mult_jobs)  
                oneCalc['__execCounterBkgrd__']+=1
		AFLOWpi.prep.__addToBlock(oneCalc,ID,'RUN', execString)
#####################################################################                

        
        return calcs

#####################################################################################################################




########################################################################################################################################################################################################################################

def construct_and_run(__submitNodeName__,oneCalc,ID,build_command='',subset_tasks=[],fault_tolerant=False,mult_jobs=True,subset_name='SUBSET',keep_file_names=False,clean_input=True):


	'''this is a check to see if we're restarting when mult_jobs==True'''
	checkBool=False
	if '__CRAWL_CHECK__' in oneCalc.keys():
            if oneCalc['__CRAWL_CHECK__']==ID:
                checkBool=True


	chain_index=1
	try:

            chain_index=oneCalc['__chain_index__']
            #		AFLOWpi.prep.__passGlobalVar('__TEMP__INDEX__COUNTER__',oneCalc['__TEMP__INDEX__COUNTER__'])
            AFLOWpi.prep.__passGlobalVar('__TEMP__INDEX__COUNTER__',chain_index)
	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
	
	chain_logname='step_%02d'%chain_index


        if checkBool==False:
            #block this prep from being run again.
            oneCalc['__CRAWL_CHECK__']=ID
            AFLOWpi.prep.__saveOneCalc(oneCalc,ID)

            outFile=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in')
            command = '''
         completeBool=True
	 if completeBool:
	    workdir = AFLOWpi.prep.__ConfigSectionMap('prep','workdir')
	    mainOneCalc = AFLOWpi.prep.__loadOneCalc(workdir,'%s')
	    AFLOWpi.run.__submitJob('%s',mainOneCalc,__submitNodeName__,forceOneJob=True)

''' % (ID,ID)
################################################################################################################
            try:
                os.mkdir(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi'))
            except:
                pass

            newConfigPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi','CONFIG.config')
            config = ConfigParser.RawConfigParser()
            config.read(oneCalc['_AFLOWPI_CONFIG_'])
            config.set('prep', 'workdir', oneCalc['_AFLOWPI_FOLDER_']) 
            if config.has_section('cluster'):
                try:

                    qsub_temp_ref = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi','CLUSTER.ref')
                    qsubSub='''cd %s\npython %s''' % (os.path.abspath(oneCalc['_AFLOWPI_FOLDER_']),os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py')))
                    qsubSub_reg = re.compile(qsubSub)
                    
                    with open(oneCalc['__qsubFileName__'],'r') as qsub_pre_trans:
                        qsub_string = qsub_pre_trans.read()

                    qsub_string = qsubSub_reg.sub('',qsub_string)

                    with open(qsub_temp_ref,'w') as qsub_post_trans:
                        qsub_post_trans.write(qsub_string)

                    config.set('cluster', 'jobreffile',qsub_temp_ref) 
                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)



            with open(newConfigPath,'w') as fileWrite:    
                config.write(fileWrite)




################################################################################################################



            calc_subset = AFLOWpi.prep.__one_test_build(oneCalc,ID,build_command,subset_name=subset_name,keep_name=keep_file_names,config=newConfigPath,clean_input=clean_input)

            AFLOWpi.prep.runAfterAllDone(calc_subset,command,faultTolerant=fault_tolerant)


            '''if we are submitting the grid calc jobs separately or one big job'''



            for task in subset_tasks:
                exec(task)
                


            for ID_new,oneCalc_new in calc_subset.iteritems():
                try:
                    calc_subset[ID_new]['__walltime_dict__']=oneCalc['__walltime_dict__']
                    AFLOWpi.prep.__saveOneCalc(oneCalc_new,ID_new)
                except Exception,e:
                    print e
                    pass
                    


            '''submit in reverse order because calcs later in the orderedDict are more likely'''
            '''to be larger cells than those at the beginning'''		
 #           invert_bool=True
            '''if we're almost at the end of the walltime don't try to submit'''
            walltime,startScript=AFLOWpi.run.__grabWalltime(oneCalc,ID)

            #keep track of time in main script as it loops in case
            #we are running serial jobs and the walltime runs out
        else:
            subset_config = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi','CONFIG.config')

            calc_subset=AFLOWpi.prep.loadlogs(subset_name,'',chain_logname,config=subset_config)
            try:
                walltime,startScript=AFLOWpi.run.__grabWalltime(oneCalc,ID)
            except:
                pass
            #exit if all calcs are done (needed for local mode)
            if AFLOWpi.prep.__checkSuccessCompletion(oneCalc,ID,faultTolerant=fault_tolerant):
                return oneCalc,ID
         
        
#        for ID_new,oneCalc_new in calc_subset.iteritems():
#            calc_subset[ID_new]['__walltime_dict__']=oneCalc['__walltime_dict__']
        if mult_jobs==True:
            oneJobBool=False
            sajO=True
        else:
            oneJobBool=True
            sajO=False


        try:
                last=len(calc_subset)
                for ID_new,oneCalc_new in reversed(calc_subset.items()):

                    last-=1 
                    if last==0:
                        #to make sure this doesn't try to run again
                        oneCalc['__execCounter__']+=1
                        oneCalc['prev'].append(ID)
                        AFLOWpi.prep.__saveOneCalc(oneCalc,ID)

                        oneJobBool=True
                        sajO=False

                    AFLOWpi.run.__submitJob(ID_new,oneCalc_new,__submitNodeName__,forceOneJob=oneJobBool,sajOverride=sajO)
                    #in the case for mult_jobs=False check to see if we're at the end time
                    if mult_jobs==False:
                        if  AFLOWpi.prep.__ConfigSectionMap('cluster','type')!='':
                            if numpy.abs(time.time()-startScript)>walltime*0.90:
                                if mult_jobs==False:
                                    AFLOWpi.run.__submitJob(ID,oneCalc,__submitNodeName__,sajOverride=True)
                                    sys.exit(0)			



                sys.exit(0)                    

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            sys.exit(0)
            
        
            
        return oneCalc,ID


import shelve
import os
import datetime
import cPickle
import logging 
import collections
import AFLOWpi
import re
import subprocess
import string
import decimal
import StringIO
import random
import shutil
import ConfigParser
import sys
import AFLOWpi.qe
import __main__
import numpy
import glob
import time
import multiprocessing
import Queue
import copy
import AFLOWpi.plot
import functools


#########################################################################################################################
#########################################################################################################################
##LOCKING WRAPPERS
#########################################################################################################################
######################################################################################################################### 

def __check_lock(func,oneCalc,ID,*args,**kwargs):
	'''
	A function that is wrapped on another function to check if the script 
	has already run through to prevent certain functions in the <ID>.py 
	from running again after a restart or a loop like scfuj.

	Arguments:
	      func (func): needed for the wrapper
	      oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
 	      ID (str): ID string for the particular calculation and step
	      *args: additional arguments 

	Keyword Arguments:
	      **kwargs: additional keyword arguments

	Returns:
	      Returns either False or the calculation dictionary and its ID hash

	'''

	if 'override_lock' in kwargs.keys():
		
		return oneCalc,ID
        try:
            if ID in oneCalc['prev']:
                return False,False
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        new_ID=ID
        oneCalc['__execCounter__']=0
        oneCalc['__qsubFileName__']=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.qsub' % new_ID) 
        try:
            d=oneCalc
            f=ID
            new_calcs={}
            output_calcs = {}
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


        a = f+'.out'

        return oneCalc,ID



def newstepWrapper(pre):
    '''
    A function that wraps another function that is to be run before the 
    a certain function runs. Its use must be in the form:
    
    @newstepwrapper(func)
    def being_wrapped(oneCalc,ID,*args,**kwargs)

    where func is the function that is to be run before and being_wrapped
    is the function being wrapped. oneCalc and ID must be the first two 
    arguments in the function being wrapped. additional arguments and 
    keyword arguments can follow.
       
    Arguments:
          pre (func): function object that is to be wrapped before another function runs

    Returns:
          the returned values of function being wrapped or the execution of the function
          is skipped entirely.
	


    '''
    def decorate(func):
        def call(*args, **kwargs):
		args=list(args)
		try:

			new_oneCalc,new_ID = pre(func, *args, **kwargs)

			if False==new_oneCalc and False==new_ID:
				return args[0],args[1]
			else:

				args[0]=new_oneCalc
		except:
			pass

		args[0],args[1] = func(*args, **kwargs)
		
		return args[0],args[1]
        return call
    return decorate


def __lock_transform(oneCalc,ID):
	"""
	To be used with AFLOWpi.prep.newstepWrapper to stop a function from being run in the 
	_<ID>.py python scripts generated by AFLOWpi to prohibit the wrapped function 
	from running if it has already run.

	Arguments:
             oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
             ID (str): ID string for the particular calculation and step

	  
	Returns:
             oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
             ID (str): ID string for the particular calculation and step
	     
	"""
#	try:
#		globals()['__TEMP__INDEX__COUNTER__']=oneCalc['__chain_index__']
#	except:
#		pass
	'''make sure it only goes through this once'''
        try:
            if ID in oneCalc['prev']:
                return oneCalc,ID
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


#	globals()['__TEMP__INDEX__COUNTER__']=oneCalc['__chain_index__']

        subdir=oneCalc['_AFLOWPI_FOLDER_']
        oneCalc['execFile']= os.path.join(subdir,'_%s.py' % ID)

        '''set prev to current ID so that we can check if we have already transformed'''
        oneCalc['prev'].append(ID)
        
        '''update status of the calc after it's transformed'''
        oneCalc['__status__']=collections.OrderedDict({"Start":False,'Complete':False,"Restart":0,"Error":'None'})

        """save it so if it's loaded by a scipt it'll have the correct vals"""
        __saveOneCalc(oneCalc,ID)

	return oneCalc,ID




#########################################################################################################################
#########################################################################################################################
##QE INPUT CLEANING
#########################################################################################################################
######################################################################################################################### 

def __resolveEqualities(inputString):
    '''
    Takes an input string that may or may not have text patterns of the form:

    ===EQN===

    where EQN is some syntactically correct simple equations like 5+3. _AFLOWpi keywords
    can be used in the reference input. An example of this could be in to generate a set 
    of calculations used to make an energy map of one layer of a layered material sliding
    across the other one might make their QE ATOMIC_POSITIONS card in their reference 
    input like so:

    ATOMIC_POSTIONS {crystal}
    _AFLOWPI_A_ ===0.00+_AFLOWPI_XSHIFT_===  ===0.00+_AFLOWPI_YSHIFT_===  0.00
    _AFLOWPI_A_ ===0.50+_AFLOWPI_XSHIFT_===  ===0.50+_AFLOWPI_YSHIFT_===  0.00
    _AFLOWPI_A_ ===0.00+_AFLOWPI_XSHIFT_===  ===0.00+_AFLOWPI_YSHIFT_===  0.00
    _AFLOWPI_A_ ===0.50+_AFLOWPI_XSHIFT_===  ===0.50+_AFLOWPI_YSHIFT_===  0.00
    _AFLOWPI_A_ 0.00                     0.00                     0.75
    _AFLOWPI_A_ 0.00                     0.50                     0.75
    _AFLOWPI_A_ 0.50                     0.00                     0.75
    _AFLOWPI_A_ 0.50                     0.50                     0.75

    With their user script containing a range of values for _AFLOWPI_XSHIFT_ and _AFLOWPI_YSHIFT_.	
    
    Arguments:
          inputString (str): a string of text

    Returns:
          a string of text that has the equations (i.e. ===EQN===) replaced with values

    '''


    #search for the pattern
    equalities = re.findall(ur'===[^,\n]*===',inputString)
    allvars={}
    #get a all the variables and their values into a dictionary

    origEquality=copy.deepcopy(equalities)
    equalDict={}


    for replacement in equalities:
        try:
            evalResult=eval(replacement[3:-3])
            replacementEval = str(evalResult)
            equalDict[replacement]=replacementEval
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

    for orig,replacement in equalDict.iteritems():
        try:
            inputString = re.sub(re.escape(orig),replacement,inputString)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


    return inputString



def __removeComments(inputString):
    '''
    Removes any text from a string that follows either a # or ! on to the end of
    that line.
    
    Arguments:
          inputString (str): input string with ! or # as comment characters


    Returns:
          the same string with the comments removed

    '''

    inputString = re.sub(r'!.*\n','\n',inputString)
    inputString = re.sub(r'#.*\n','\n',inputString)
    return inputString

def remove_blank_lines(inp_str):
    '''
    Removes whitespace lines of text
    
    Arguments:
          inp_str (str): input string of text

    Returns:
          the same string with blank lines removed

    '''
	
    out_str=''
    for line in inp_str.split('\n'):
        if len(line.strip())!=0:
            out_str+=line+'\n'

    return out_str



def __cleanInputStringSCF(inputString):
	'''
	Sanitizes QE input files by removing whitespace, commented out text, or any other
	unnecessary characters.
    
	Arguments:
	      inputString (str): a string of a QE input file

	Returns:
	      a string of the sanitized QE input file

	'''

        logging.debug('Entering __cleanInputString')

        removeBlank = []
        inputString = AFLOWpi.prep.__removeComments(inputString)

	inputString = AFLOWpi.prep.remove_blank_lines(inputString)
        splitInput = inputString.split('\n')         

        for lineNum in range(len(splitInput)): 
           if len(re.findall(r'\s*&control\s*',splitInput[lineNum]))==0:
               splitInput[lineNum]=''
           else:
                break
        splitLines = []
        for lines in splitInput:
            try:
                oneLine = lines.split('!')[0]            
            except:
                oneLine=lines
            
            if len(oneLine)!=0:
                splitLines.append(oneLine)

        inputString = '\n'.join(splitLines)
        
        atomicSpecRegex = re.compile(r'ATOMIC_SPECIES\s*\n((?:(?:\s*[A-Za-z0-9]+)\s+(?:[0-9.]+)\s+(?:.+)\n*)+)(?=(?:[A-Z_\s]+\n)|)',(re.MULTILINE))

	speciesSortList = []
	try:
            speciesArr = atomicSpecRegex.findall(inputString)
            speciesArr = speciesArr[-1]
            splitSpecies = speciesArr.split('\n')
	except:
            splitSpecies=''
                                                                                                                                                                         
        try:
            for item in splitSpecies:
		if item!='' and item[0]!='!':

                    speciesSplit = ' '.join([word.strip()+' ' for word in item.split()])
                    
                    if speciesSplit not in speciesSortList:
                        speciesSortList.append(speciesSplit)
#            speciesSortList = sorted(speciesSortList)
        
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        
        testDictString='ATOMIC_SPECIES\n'        

	for speciesSplit in speciesSortList:
		testDictString+=speciesSplit+'\n'
        logging.debug(testDictString)
        try:
            cleanedInputString = atomicSpecRegex.sub(r'%s\n' % testDictString,inputString)

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


        try:
           cleanedInputString = re.sub(r'ntyp\s*=\s*\d+','ntyp = %s' % len(speciesSortList),cleanedInputString)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

#        atomicPosRegex = re.compile(r"(^ *[A-Z][a-z]*)[\s]+[\d\.-]+\s+[\d\.-]+\s+[\d\.-]+", re.MULTILINE)        
        atomicPosRegex =  AFLOWpi.qe.regex.atomic_positions(cleanedInputString,'content','regex')

        pos,flag = AFLOWpi.retr.detachPosFlags(AFLOWpi.qe.regex.atomic_positions(cleanedInputString))
        try:
            posList =  AFLOWpi.qe.regex.atomic_positions(cleanedInputString)
            posList,flag = AFLOWpi.retr.detachPosFlags(posList)
            posList=posList.split('\n')

            
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        try:
            posList.remove('')
        except:
            pass
        try:
            posList.remove('\r')
        except:
            pass
        try:
            posList.remove('\n')
        except:
            pass

        for entry in range(len(posList)):
            try:
                posList[entry]=' '.join([x.rstrip() for x in  posList[entry].split()])
            except Exception,e:
                print e

        posListNew = []
        
        posListNewString = '\n'.join(posList)
        posListNewString = AFLOWpi.retr.attachPosFlags(posListNewString,flag)
        nat =  len(posList)

        

        cleanedInputString = re.sub(r'nat\s*=\s*\d+','nat = %s' % nat,cleanedInputString)


        cellParamRegex =  AFLOWpi.qe.regex.atomic_positions(cleanedInputString,'modifier','regex')
        if len(cellParamRegex.findall(inputString)):
            paramUnits=cellParamRegex.findall(inputString)[-1].strip()
            paramUnits = paramUnits.lower()
            

        splitInput = AFLOWpi.retr.__splitInput(cleanedInputString)
        splitInput['ATOMIC_POSITIONS']['__content__']=posListNewString 
        splitInput['ATOMIC_POSITIONS']['__modifier__']=paramUnits
        cleanedInputString  = AFLOWpi.retr.__joinInput(splitInput)
        

        cleanedInputString=AFLOWpi.prep.__transformInput(cleanedInputString)

        logging.debug('Exiting __cleanInputString')


        return cleanedInputString







def __parseRef(refFile,dictFlag=False):
	'''
	Parses the input reference file and extract all of the information
	it then organizes it in a standard way and is used to create a 
	unique hash checksum for the fortran namelist input. This is so spacing
	and order of namelist input data will not affect the checksum of each
	unique calculation input file

	Arguments:
 	      refFile (str): a string of the input reference file

	Keyword Arguments:
	      dictFlag (bool): whether you want to have a string of the cleaned up text or a 
   	                 tokenized dictionary of it.
	

	Returns: 
    	      Dictionary or String
	
        '''


        refFile = AFLOWpi.prep.__removeComments(refFile)

	refFile = AFLOWpi.prep.remove_blank_lines(refFile)
        refFileSplit = AFLOWpi.retr.__splitInput(refFile)

        refFileSplit = AFLOWpi.retr.__orderSplitInput(refFileSplit)
        del refFileSplit['&control']
        refFile = AFLOWpi.retr.__joinInput(refFileSplit)
        refFileSplit = AFLOWpi.retr.__splitInput(refFile)
#	cellParamListRegex = re.compile(r'(ATOMIC_POSITIONS\s*.+\n)')
#        icPosRegex =	AFLOWpi.qevars.atomic_positions

	######################
	####ATOMIC SPECIES####
	######################
        try:
            try:
                speciesArr=refFileSplit['ATOMIC_SPECIES']['__content__']
                splitSpecies = speciesArr.split('\n')
            except ValueError:
                pass
            speciesSortList=[]
            for item in splitSpecies:
                if item!='':
                    speciesSplit = ' '.join(item.split())
                    speciesSortList.append(speciesSplit)
                    speciesSortList = sorted(speciesSortList)
                    speciesSplit = '\n'.join(speciesSortList)
                    refFileSplit['ATOMIC_SPECIES']['__content__'] = speciesSplit
        
               
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

	######################
	###ATOMIC POSITIONS###
	######################
        try:
            splitAtomicPos= [x for x in refFileSplit['ATOMIC_POSITIONS']['__content__'].split('\n') if len(x.strip())]
            splitAtomicPos,flags=AFLOWpi.retr.detachPosFlags('\n'.join(splitAtomicPos))
            atomPosSortList=[]
            for item in splitAtomicPos:
                    if item!='':
                            atomSplit = item.split()
                            formattedAtomPos = '%3.3s' % atomSplit[0].strip()+' '+' '.join(['%16.14f' % (decimal.Decimal(str(atomSplit[i]))) for i in range(1,len(atomSplit))])

                            
                            atomPosSortList.append(formattedAtomPos)

            flagsSplit= [x.strip() for x in flags.split('\n')]
#            formatFlags = ' '.join([' '.join([x.strip() for x in flagSplit[i].split()]) for i in range(len(flagSplit))])

            atomPosSortList = sorted(atomPosSortList)
            positionSplit = '\n'.join(atomPosSortList)
            refFileSplit['ATOMIC_POSITIONS']['__content__'] = positionSplit
        except:
            pass
        refFile = AFLOWpi.retr.__joinInput(refFileSplit)

	###################
	####CELL PARAMS####
	###################

        refFileSplit = AFLOWpi.retr.__splitInput(refFile)
        for item in refFileSplit['&system'].keys():
            if len(re.findall('celldm',item))!=0:
                celldm = refFileSplit['&system'][item]
                refFileSplit['&system'][item] = '%14.12f' % numpy.around(numpy.float(celldm),decimals=5)

        kpointString = refFileSplit['K_POINTS']['__content__']
        kpointMod = refFileSplit['K_POINTS']['__modifier__']


        kpointSplit =  [x for x in kpointString.split('\n') if len(x)!=0]

        kpointSplitFixed=[]
        for x in kpointSplit:
            splitX=[]
            if len(x)!=0:
                splitX=x.split()
                for y in range(len(splitX)):
                    try:
                        splitX[y]='%1.16f' % float(splitX[y])
                    except:
                        pass
                kpointSplitFixed.append(' '.join(splitX))

        kpointString = '\n'.join(sorted(kpointSplitFixed))

        refFileSplit['K_POINTS']['__content__'] = kpointString

        
        
        refFile = AFLOWpi.retr.__joinInput(refFileSplit)


        '''sorting each of the namelists'''
        refFileSplit = AFLOWpi.retr.__splitInput(refFile)
        for namelist,items in refFileSplit.iteritems():
            
            if '&' in namelist:
                refFileSplit[namelist]=collections.OrderedDict(sorted(items.items()))
        

        refFile = AFLOWpi.retr.__joinInput(refFileSplit)



        return refFile

        #######################################################################################################
        ####TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED ###
        #######################################################################################################

	###################
	###OCCUPATIONS###
	###################
	OCCNLRegex      = re.compile(r'(OCCUPATIONS\s*\n)',re.MULTILINE)
	OCCRegex     = re.compile(r'(?:OCCUPATIONS\s*\n)((?:[A-Za-z0-9|\s|.]+\n*)+)(?=(?:[A-Z|_|\s]+\n)|)',re.MULTILINE)
	OCCSortList = []
	try:
		try:

			inputStr+=OCCNLRegex.findall(refFile)[-1].strip()+'\n'
		except ValueError:
			pass
	
		try:
			OCCArr = OCCRegex.findall(refFile)[-1]

			splitOCC = OCCArr.split('\n')
		except ValueError:
			pass
		for item in splitOCC:
			if item!='':
				atomSplit = item.split()
				formattedOCC = ' '.join(['%g' % (decimal.Decimal(str(atomSplit[i]))) for i in range(len(atomSplit))])+'\n'
				OCCSortList.append(formattedOCC)
		OCCSortList = sorted(OCCSortList)
		testDictString=''
		for formattedOCC in OCCSortList:
			inputStr+=formattedOCC
			testDictString+=formattedOCC
		inDict['_AFLOWPI_CELLPARAM_']=testDictString
			
	except:
		pass


	###################
	####CONSTRAINTS####
	###################
        CONSTNLRegex      = re.compile(r'(CONSTRAINTS\s*\n)',re.MULTILINE)
	CONSTRegex     = re.compile(r"(?:CONSTRAINTS\s*\n)((?:[a-z0-9\s.'_]+\n*)+)(?=(?:[A-Z_\s]+\n)|)",re.MULTILINE)
	CONSTSortList = []
	try:
		try:

			inputStr+=CONSTNLRegex.findall(refFile)[-1].strip()+'\n'
		except ValueError:
			pass
	
		try:
			CONSTArr = CONSTRegex.findall(refFile)[-1]

			splitCONST = CONSTArr.split('\n')
		except ValueError:
			pass
		for item in splitCONST:
			if item!='':
				
				atomSplit = item.split()
				
				formattedCONST = atomSplit[0]+' '+' '.join(['%g' % (decimal.Decimal(str(atomSplit[i]))) for i in range(1,len(atomSplit))])+'\n'
				CONSTSortList.append(formattedCONST)
		CONSTSortList = sorted(CONSTSortList)
		testDictString=''
		for formattedCONST in CONSTSortList:
			inputStr+=formattedCONST
			testDictString+=formattedCONST
		inDict['_AFLOWPI_CONSTRAINTS_']=testDictString
	except Exception,e:
		pass


	
	###################
	###ATOMIC FORCES###
	###################

	atomicForNLRegex      = re.compile(r'(ATOMIC_FORCES\s*\n)',re.MULTILINE)
	atomicForRegex     = re.compile(r'(?:ATOMIC_FORCES\s*\n)((?:[A-Za-z0-9]+\s+(?:(?:[0-9|.]+)\s*){3}\n*)+)(?=(?:[A-Z|_|\s]+\n)|)',re.MULTILINE)
	atomForSortList = []
	try:
		try:

			inputStr+=atomicForNLRegex.findall(refFile)[-1].strip()+'\n'
		except ValueError:
			pass
	
		try:
			atomicForArr = atomicForRegex.findall(refFile)[-1]

			splitAtomicFor = atomicForArr.split('\n')
		except ValueError:
			pass
		for item in splitAtomicFor:
			if item!='':
				atomSplit = item.split()
				formattedAtomFor = atomSplit[0]+' '+' '.join(['%g' % (decimal.Decimal(str(atomSplit[i]))) for i in range(1,len(atomSplit))])+'\n'
				atomForSortList.append(formattedAtomFor)
		atomForSortList = sorted(atomForSortList)
		testDictString=''
		for formattedAtomFor in atomForSortList:
			inputStr+=formattedAtomFor
			testDictString+=formattedAtomFor
		inDict['_AFLOWPI_INPUTFORCES_']=testDictString
	except:
		pass





        logging.debug('exiting __parseRef')	
	if not dictFlag:
		return inputStr
	else:
		return inDict

        #######################################################################################################
        ####TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED TO BE COMPLETED ###
        #######################################################################################################


def __transformParamsInput(inputString):
    '''
    Transforms the various styles of defining the lattice vectors into celldm form so AFLOWpi
    has a standard format to work with.

	
    Arguments:
          inputString (str): a string of a QE input

    Returns:
          a string of a QE input that has its lattice vectors in celldm style

    '''

    inputDict = AFLOWpi.retr.__splitInput(inputString)
    BohrToAngstrom = 0.529177249

    #    and not 'celldm(1)' in inputDict['&system'].keys() and not 'A' in inputDict['&system'].keys():
    CELL_PARAM_FLAG=False


    if 'CELL_PARAMETERS' in inputDict.keys():
        cellParamMatrix = AFLOWpi.retr.__cellStringToMatrix(inputDict['CELL_PARAMETERS']['__content__'])
        cellParamRegex = AFLOWpi.qe.regex.cell_parameters('','content','regex')

        if len(cellParamRegex.findall(inputString)):

            paramUnits = inputDict['CELL_PARAMETERS']['__modifier__'].replace('{','').replace('}','')
            paramUnits = paramUnits.lower()

            if paramUnits=='angstrom':
                cellParamMatrix*=1/0.529177249

            else:
                if 'a' in inputDict['&system'].keys():
                    mult = float(inputDict['&system']['a'])
                elif paramUnits=='angstrom':
                    mult=1/0.529177249
                elif 'celldm(1)' in inputDict['&system'].keys():
                    mult = float(inputDict['&system']['celldm(1)'])
                else: 
                    mult=1.0
                cellParamMatrix*=mult

            ibravDict = AFLOWpi.retr.__free2celldm(cellParamMatrix)
	    ibravDict_prime = AFLOWpi.retr.__getCelldm2freeDict(ibravDict)
	    new_prim_mat = AFLOWpi.retr.celldm2free(**ibravDict_prime) 
	    new_prim_mat = AFLOWpi.retr.__cellStringToMatrix(new_prim_mat)
	    orig=  AFLOWpi.retr.__prim2ConvVec(cellParamMatrix)        
	    new =  AFLOWpi.retr.__prim2ConvVec(new_prim_mat)        
	    slam= numpy.linalg.inv(new.astype(numpy.float32)).dot(orig)
            ibravDict = AFLOWpi.retr.__free2celldm(cellParamMatrix)
	    
	    
            for param,value in ibravDict.iteritems():
                if param.strip()=='celldm(1)':
                    scaled = str(float(value)/BohrToAngstrom)
                    if 'a' in inputDict['&system'].keys():
                        value = str(float(value)/BohrToAngstrom)
                    inputDict['&system'][param]=str(float(value))
                else:
                    inputDict['&system'][param]=value

            try:
                del inputDict['CELL_PARAMETERS']
            except Exception,e:
                print e

            inputString = AFLOWpi.retr.__joinInput(inputDict)

            return inputString    


    elif 'A' in [x.upper() for x in inputDict['&system'].keys()]:
        paramDict={}
        A=float(inputDict['&system']['a'])
        paramDict['ibrav']=int(inputDict['&system']['ibrav'])
        inputDictCopy = copy.deepcopy(inputDict)

        for param in ['A','B','C','cosAB','cosAC','cosBC']:
                if param.lower() in inputDict['&system'].keys():
                    del inputDictCopy['&system'][param.lower()]
                    
        tempDict=collections.OrderedDict()
        for param in ['A','B','C','COSAB','COSAC','COSBC']:

                if param.lower() in inputDict['&system'].keys():
                    paramFloat = float(inputDict['&system'][param.lower()])

		    param=param.upper()
                    if param=='COSAB':
                        paramName='celldm(6)'
                    if param=='COSAC':
                        paramName='celldm(5)'
                    if param=='COSBC':
                        paramName='celldm(4)'
                    if param=='A':
                        paramName='celldm(1)'
                        paramFloat=numpy.around((paramFloat/0.529177249),decimals=5)
                    if param=='B':
                        paramName='celldm(2)'
                        paramFloat/=A
                    if param=='C':
                        paramName='celldm(3)'
                        paramFloat/=A
                    

                    tempDict.update({paramName:paramFloat})



        for k,v in tempDict.iteritems():
            inputDictCopy['&system'][k]=v




        inputString = AFLOWpi.retr.__joinInput(inputDictCopy)


        return inputString    


    elif 'celldm(1)' in inputDict['&system'].keys():
        inputDict = AFLOWpi.retr.__splitInput(inputString)        
        return AFLOWpi.retr.__joinInput(inputDict)
#         BohrToAngstrom = 0.529177249
#         paramList = ['celldm(1)','celldm(2)','celldm(3)','celldm(4)','celldm(5)','celldm(6)']
#         paramDict={}
#         paramDict['ibrav']=int(inputDict['&system']['ibrav'])
#         for param in range(len(paramList)):
#             try:
#                 paramFloat = float(inputDict['&system'][paramList[param]])
#                 print param
#                 if param<1:
#                     paramFloat*=BohrToAngstrom
#                 print paramFloat
#                 paramName = 'celldm'+str(param+1)
#                 paramDict[paramName]=paramFloat
#                 del inputDict['&system'][paramList[param]]
#             except:
#                 pass
#         paramDict['returnString']=False

#         cellParamMatrix = AFLOWpi.retr.celldm2free(**paramDict)
#         print paramDict
#         print cellParamMatrix


# #        inputDict['&system']['ibrav']='0'
#         string=''
#         for entry in cellParamMatrix.getA():
#             string+= ' '.join(['%20.14f' % (Decimal(str(entry[i]))) for i in range(len(entry))])+'\n'         

#         string+='\n'
#    inputDict['CELL_PARAMETERS (angstrom)']=string

            


def __transformPositionsInput(inputString):
    '''
    Transforms the various styles of defining the atomic positions into crystal 
    fractional coordinates form so AFLOWpi has a standard format to work with.

	
    Arguments:
          inputString (str): a string of a QE input

    Returns:
          inputString (str): a string of a QE input that has its atomic positions in crystal fractional coordinates

    '''

    cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(inputString)


    splitInput =  AFLOWpi.retr.__splitInput(inputString)
    try:

        labels =  AFLOWpi.retr.__getPosLabels(inputString)
        positionMatrix = AFLOWpi.retr.__getPositions(inputString)
        pos,flags = AFLOWpi.retr.detachPosFlags(AFLOWpi.qe.regex.atomic_positions(inputString))

        paramUnits = splitInput['ATOMIC_POSITIONS']['__modifier__']
        cellParamMatrix/=float(splitInput['&system']['celldm(1)'])


        if paramUnits=='{angstrom}':
            positionMatrix/=0.52917721092
	    positionMatrix/=float(splitInput['&system']['celldm(1)'])
            symMatrix = AFLOWpi.retr.__convertFractional(positionMatrix,cellParamMatrix)
        elif paramUnits=='{bohr}':
            positionMatrix/=float(splitInput['&system']['celldm(1)'])
            symMatrix = AFLOWpi.retr.__convertFractional(positionMatrix,cellParamMatrix)
        elif paramUnits=='{crystal}' or paramUnits=='':

            symMatrix = positionMatrix
        elif paramUnits=='alat':
            pass
        else:
            symMatrix = positionMatrix

        outputString = ''

        for entry in range(len(symMatrix.getA())):

            posLineStr = ' '.join(['%20.14f' % (decimal.Decimal(str(numpy.around(i,9)))) for i in symMatrix.getA()[entry]])+'\n'
            outputString+='%4s %8s' % (labels[entry],posLineStr)
        
        outputString = AFLOWpi.retr.attachPosFlags(outputString,flags)

        splitInput['ATOMIC_POSITIONS' ]['__content__'] = outputString
        splitInput['ATOMIC_POSITIONS' ]['__modifier__'] = '{crystal}'

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e) 

        return inputString

    returnString =   AFLOWpi.retr.__joinInput(splitInput)


    return returnString


def __transformInput(inputString):
    '''
    Standardizes the input format of the QE input files so AFLOWpi only has to account
    for one standard input type for atomic positions and lattice parameters.
    
    Arguments:
          inputString (str): a string of a QE input file

    Returns:
          inputString (str): a string of a QE input file

    '''

    try:
	    positionMatrix = AFLOWpi.retr.__getPositions(inputString)
	    labels = AFLOWpi.retr.__getPosLabels(inputString)
	    pos,flag = AFLOWpi.retr.detachPosFlags(AFLOWpi.qe.regex.atomic_positions(inputString))
    except Exception,e:
	    print e
    
    try:
	    if AFLOWpi.prep.__ConfigSectionMap('prep','convention').upper()=='AFLOW':
		    pass
#		    inputString=AFLOWpi.prep.__aflow2pw(inputString)    
    except Exception,e:
	    print e
    

    inputString = __transformParamsInput(inputString)
    inputString = __transformPositionsInput(inputString) 




    return inputString


#########################################################################################################################
#########################################################################################################################
##QE ATOMIC_SPECIES LISTS
#########################################################################################################################
######################################################################################################################### 



import os

def __getAMass(atom):
	"""
	Get the atomic mass for the specific atomic species in input
	
	Arguments:
	      atom (str): Atomic Species you want the mass of

	Returns:
	      Mass of the atom's species

	"""

	logging.debug('Entering getAMass')	
	mass = {
        
	"H": (1.008, 'Hydrogen'),
	"He": (4.003, 'Helium'),
	"Li": (6.938, 'Lithium'),
	"Be": (9.012, 'Beryllium'),
	"B": (10.806, 'Boron'),
	"C": (12.010, 'Carbon'),
	"N": (14.006, 'Nitrogen'),
	"O": (15.999, 'Oxygen'),
	"F": (18.998, 'Fluorine'),
	"Ne": (20.180, 'Neon'),
	"Na": (22.990, 'Sodium'),
	"Mg": (24.304, 'Magnesium'),
	"Al": (26.982, 'Aluminium'),
	"Si": (28.084, 'Silicon'),
	"P": (30.974, 'Phosphorus'),
	"S": (32.059, 'Sulfur'),
	"Cl": (35.446, 'Chlorine'),
	"Ar": (39.948, 'Argon'),
	"K": (39.098, 'Potassium'),
	"Ca": (40.078, 'Calcium'),
	"Sc": (44.955, 'Scandium'),
	"Ti": (47.867, 'Titanium'),
	"V": (50.941, 'Vanadium'),
	"Cr": (51.996, 'Chromium'),
	"Mn": (54.938, 'Manganese'),
	"Fe": (55.845, 'Iron'),
	"Co": (58.933, 'Cobalt'),
	"Ni": (58.693, 'Nickel'),
	"Cu": (63.546, 'Copper'),
	"Zn": (65.380, 'Zinc'),
	"Ga": (69.723, 'Gallium'),
	"Ge": (72.630, 'Germanium'),
	"As": (74.922, 'Arsenic'),
	"As": (74.922, 'Arsenic'),
	"Se": (78.971, 'Selenium'),
	"Br": (79.901, 'Bromine'),
	"Kr": (83.798, 'Krypton'),
	"Rb": (85.468, 'Rubidium'),
	"Sr": (87.620, 'Strontium'),
	"Y": (88.906, 'Yttrium'),
	"Zr": (91.224, 'Zirconium'),
	"Nb": (92.906, 'Niobium'),
	"Mo": (95.950, 'Molybdenum'),
        "Tc": (98.0, 'Tecnetium'),
	"Ru": (101.070,' Ruthenium'),
	"Rh": (102.906, 'Rhodium'),
	"Pd": (106.420, 'Palladium'),
	"Ag": (107.868, 'Silver'),
	"Cd": (112.414, 'Cadmium'),
	"In": (114.818, 'Indium'),
	"Sn": (118.710, 'Tin'),
	"Sb": (121.760, 'Antimony'),
	"Te": (127.600, 'Tellurium'),
	"I": (126.904, 'Iodine'),
	"Xe": (131.293, 'Xenon'),
	"Cs": (132.905, 'Cesium'),
	"Ba": (137.327, 'Barium'),
	"La": (138.905, 'Lanthanum'),
	"Ce": (140.116, 'Cerium'),
	"Pr": (140.908, 'Praseodymium'),
	"Nd": (144.242, 'Neodymium'),
	"Sm": (150.360, 'Samarium'),
	"Eu": (151.964, 'Europium'),
	"Gd": (157.250, 'Gadolinium'),
	"Tb": (158.925, 'Terbium'),
	"Dy": (162.500, 'Dysprosium'),
	"Ho": (164.930, 'Holmium'),
	"Er": (167.259, 'Erbium'),
	"Tm": (168.934, 'Thulium'),
	"Yb": (173.054, 'Ytterbium'),
	"Lu": (174.967, 'Lutetium'),
	"Hf": (178.490, 'Hafnium'),
	"Ta": (180.947, 'Tantalum'),
	"W": (183.840, 'Tungsten'),
	"Re": (186.207, 'Rhenium'),
	"Os": (190.230, 'Osmium'),
	"Ir": (192.217, 'Iridium'),
	"Pt": (195.084, 'Platinum'),
	"Au": (196.966, 'Gold'),
	"Hg": (200.592, 'Mercury'),
	"Tl": (204.382, 'Thallium'),
	"Pb": (207.200, 'Lead'),
	"Bi": (208.980, 'Bismuth'),
        "Po": (209, 'Polonium'),
        "At": (210, 'Astatine'),
	"Th": (232.038, 'Thorium'),
	"Pa": (231.036, 'Protactinium'),
	"U": (238.029, 'Uranium'),
	"Am": (243.000, 'Americium'),
        "!": (0.000, 'VACANT'),}


	logging.debug('Exiting getAMass')
        atomMass=0
        try:
            atomMass = mass[atom][0]

            return atomMass
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            return 0

def __getPseudofilename(atom,pseudodir=os.path.join(os.curdir,'PSEUDOs')):
	"""
	Gets the pseudopotential filename for the specific atomic species in input.
	It will check in the AFLOWpi config file used when itiating the session for 
	a pseudodir path.

	Species names must be in the form of the element's symbol possibly followed by integers
	(i.e) "O1","Co6", or "Sn165" (Note that QE needs to be modified to accept species labels
	longer than 3 characters. Instructions on how are in the engine_mods directory.)
	
	Arguments:
	      atom (str): a string designating the atomic species you want the pseudofile name for
	
	Keyword Arguments:
	      pseudodir (str): the path of the directory containing pseudofiles

	Returns:
	      pseudofilename (str): name of the pseudopotential file in for 
	                            that species in the pseudodir specified
	      
	"""

	logging.debug('Entering getPseudofilename')
	if pseudodir == None:
	        pseudodir= __ConfigSectionMap('prep','pseudodir')
                
                if os.path.isabs(pseudodir) == False:
                    configFileLocation = __getConfigFile()



        if atom != '!':
            atom0=''
            atom1=''
            atom0 = atom[0]
            try:
                atom1 = atom[1]
            except:
                atom1 = ''
            species='['+atom0.upper()+atom0.lower()+']'+'['+atom1.upper()+atom1.lower()+']'
            regexFromConfig = __ConfigSectionMap('prep','pseudoregex')
            if regexFromConfig=='':
                regexFromConfig="re.compile(r'^'+species+'[._-]', re.I)"
            rex1 = eval(regexFromConfig)

            for l in os.listdir(pseudodir):
                    if rex1.search(l):
                            pseudofilename = l
                            logging.debug('Exiting getPseudofilename')
                            return pseudofilename
            logging.error('Missing Pseudopotential File')
            logging.debug('Exiting getPseudofilename')
            return None
        else:
            return ''







#########################################################################################################################
#########################################################################################################################
## step preparation
#########################################################################################################################
######################################################################################################################### 

            
def writeToScript(executable,calcs,from_step=0):
    """
    Generates calls on several functions to set up everything that is needed for a new step in the workflow.
    The mechanics of the _ID.py are written to it here. 

    Arguments:
	  calcs (dict): dictionary of dictionaries of calculations
          executable (str): <DEFUNCT OPTION: HERE FOR LEGACY SUPPORT> 
	  *args: <DEFUNCT OPTION: HERE FOR LEGACY SUPPORT>

    Keyword Arguments:
	  **kwargs: <DEFUNCT OPTION: HERE FOR LEGACY SUPPORT>

    Returns:
          A set of calculations for a new step in the workflow

    """

    new_calcs = collections.OrderedDict()
#    if '__TEMP__INDEX__COUNTER__' not in globals().keys():
#        globals()['__TEMP__INDEX__COUNTER__']=1
#    else:
#        globals()['__TEMP__INDEX__COUNTER__']+=1        

    extension=''

    for ID,oneCalc in calcs.iteritems():
        new_oneCalc,new_ID = __writeToScript(executable,oneCalc,ID,from_step=from_step)
    
        new_oneCalc['__status__']=collections.OrderedDict({"Start":False,'Complete':False,"Restart":0,"Error":'None'})    
	
	
        if AFLOWpi.prep.__ConfigSectionMap("cluster","type") != '':              
            AFLOWpi.run.__qsubGen(new_oneCalc,new_ID)     
            new_oneCalc['__qsubFileName__']=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s%s.qsub' % (new_ID,extension) )
#        __saveOneCalc(new_oneCalc,new_ID)
        new_calcs[new_ID]=new_oneCalc
	step_index=new_oneCalc['__chain_index__']
    
    updatelogs(new_calcs,'step_%02d'%step_index,runlocal=True)
    return new_calcs




def __fillTemplate(oneCalc,ID):
    """
    Fills in the blocks of each calculation's _ID.py with mechanism for starting the next step in the chain
    

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	  ID (str): ID string for the particular calculation and step

    Returns:
          None

    """

    logging.debug('entering __fillTemplate')
#    try:
#	    null=globals()["__TEMP__INDEX__COUNTER__"]
#    except:
#	    globals()["__TEMP__INDEX__COUNTER__"]=1
    try:
        logfile = os.path.abspath(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'AFLOWpi','LOG.log')) 
        
        if __findInBlock(oneCalc,ID,'LOGGING',"import logging") == False:
            __addToBlock(oneCalc,ID,'LOGGING',"import logging")
            __addToBlock(oneCalc,ID,'LOGGING',"logfile='%s'" % logfile)
            __addToBlock(oneCalc,ID,'LOGGING',"logging.basicConfig(filename=logfile,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',level=AFLOWpi.prep.__getLoglevel())\n\n")

        if __findInBlock(oneCalc,ID,'ID',"ID='%s'" % ID)==False:
            __addToBlock(oneCalc,ID,'ID',"ID='%s'" % ID)

        if __findInBlock(oneCalc,ID,'IMPORT',"import AFLOWpi") == False:
            __addToBlock(oneCalc,ID,'IMPORT',"import AFLOWpi")
            __addToBlock(oneCalc,ID,'IMPORT',"import time")
#            __addToBlock(oneCalc,ID,'IMPORT',"if '__GLOBAL_CLOCK__' not in globals().keys():")

            __addToBlock(oneCalc,ID,'IMPORT','__CLOCK__=time.time()')
	    __addToBlock(oneCalc,ID,'IMPORT',"import os")
	    __addToBlock(oneCalc,ID,'IMPORT',"os.chdir('%s')"%oneCalc['_AFLOWPI_FOLDER_'])
	    
#            __addToBlock(oneCalc,ID,'IMPORT','   globals()["__GLOBAL_CLOCK__"]=__GLOBAL_CLOCK__')
#            __addToBlock(oneCalc,ID,'IMPORT','else:')
#            __addToBlock(oneCalc,ID,'IMPORT','   __GLOBAL_CLOCK__ = globals()["__GLOBAL_CLOCK__"]')
            clusterType  = __ConfigSectionMap('cluster','type').lower()
            if clusterType!='':

#                __addToBlock(oneCalc,ID,'IMPORT',"from multiprocessing import Process")
#################################################################################################################
#################################################################################################################
#### NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED  ####
#################################################################################################################
#################################################################################################################
#                __addToBlock(oneCalc,ID,'RESTART',"AFLOWpi.prep.__to_local_scratch(oneCalc,ID)")
                __addToBlock(oneCalc,ID,'IMPORT',"from multiprocessing import os,Process")


                __addToBlock(oneCalc,ID,'RESTART',"oneCalc = AFLOWpi.run.__setStartTime(oneCalc,ID)")
                __addToBlock(oneCalc,ID,'RESTART',"__restartDaemon__ = Process(target=AFLOWpi.run.__restartScript, args=(oneCalc,ID,os.getpid(),))")
                __addToBlock(oneCalc,ID,'RESTART',"__restartDaemon__.daemon=True")
#                __addToBlock(oneCalc,ID,'RESTART',"__restartDaemon__.start()")
                __addToBlock(oneCalc,ID,'RESTART',"import signal")
                __addToBlock(oneCalc,ID,'RESTART',"signal.signal(signal.SIGUSR1,AFLOWpi.run.__exitClean)")       
                __addToBlock(oneCalc,ID,'RESTART',"signal.signal(signal.SIGTERM,AFLOWpi.run.__recordDeath)")       





#################################################################################################################
#################################################################################################################
#### NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED ##NOT FINISHED  ####
#################################################################################################################
#################################################################################################################
        '''switch over to local config file in AFLOWpi folder of calc tree'''
        configFile = os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'AFLOWpi','CONFIG.config')

        if __findInBlock(oneCalc,ID,'CONFIGFILE',"configFile") == False:
            __addToBlock(oneCalc,ID,'CONFIGFILE',"configFile='%s'" % configFile)
            __addToBlock(oneCalc,ID,'CONFIGFILE',"AFLOWpi.prep.__forceGlobalConfigFile(configFile)")

        try:
            if '__submitNodeName__' not in globals().keys():

                global __submitNodeName__
                __submitNodeName__ = socket.gethostname()

            if __findInBlock(oneCalc,ID,'CONFIGFILE',"__submitNodeName__ = '%s'" % __submitNodeName__)==False:
                __addToBlock(oneCalc,ID,'CONFIGFILE',"__submitNodeName__ = '%s'" % __submitNodeName__)
                __addToBlock(oneCalc,ID,'CONFIGFILE',"AFLOWpi.prep.__forceSubmitNodeIP(__submitNodeName__)")


            index=1
            try:
#                index =int(re.findall('step_(\d*)_*',ID)[-1])
		index=int(oneCalc['__chain_index__'])
            except:
                pass

#            if index==1:
#                fileName = os.path.abspath(os.path.join((os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])),'AFLOWpi','calclogs','calclog.log'))
#            else:
	    fileName = os.path.abspath(os.path.join((os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])),'AFLOWpi','calclogs','calclog_%.02d.log' % (index)))
            if __findInBlock(oneCalc,ID,'POSTPROCESSING',"AFLOWpi.prep.__updatelogs(ID,oneCalc,'%s')" % fileName)==False:
                #make sure the restart mode is reset before exiting to the next step
                __addToBlock(oneCalc,ID,'POSTPROCESSING',"""oneCalc,ID=AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','restart_mode','"from_scratch"')""")
                __addToBlock(oneCalc,ID,'POSTPROCESSING',"""AFLOWpi.prep.__updatelogs(ID,oneCalc,'%s')""" % fileName)





#            __checkSuccessCompletion({ID:oneCalc})


        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    
    logging.debug('exiting __fillTemplate')



def __getNextOneCalcVarName(oneCalc,ID):
    '''
    ##############
    ###OBSOLETE###
    ##############
    Gives increased the ID of one calculation an increase of one in its postfix 

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	  ID (str): ID string for the particular calculation and step

    Returns:
          the identical calculation object that was inputted and its ID that has had its ID modified

    '''


    logging.debug('entering __getNextOneCalcVarName')
    filename = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.py' % ID)


    try:
        with open(filename,'r') as scriptFile:
            scriptFileString = scriptFile.read()
        varList = set(re.findall(r'oneCalc',scriptFileString))    

    except:
        return 'ID','oneCalc'

    logging.debug('exiting __getNextOneCalcVarName')

    nextIDName = 'ID_%s' % (len(varList)+1)
    nextCalcName = 'oneCalc_%s' % (len(varList)+1)
    return nextIDName,nextCalcName



def __writeTemplate(finp):
    '''
    Writes the blocks of the skeleton _ID.py to file for one calculation.

    Arguments:
          finp (str): filepath for skeleton _ID.py file to go
	 
    Returns:
          None

    '''

    try:
        with open(finp,'w') as fileName:
            fileName.write('''
#####AFLOWpi EXECUTABLE#######

#IMPORT_BLOCK

#END_IMPORT_BLOCK
#CONFIGFILE_BLOCK

#END_CONFIGFILE_BLOCK
#LOGGING_BLOCK

#END_LOGGING_BLOCK
#ONECALC_BLOCK

#END_ONECALC_BLOCK
#LOADCALC_BLOCK

#END_LOADCALC_BLOCK
#ID_BLOCK

#END_ID_BLOCK
#RESTART_BLOCK

#END_RESTART_BLOCK
#PREPROCESSING_BLOCK

#END_PREPROCESSING_BLOCK

#LOCK_BLOCK

#END_LOCK_BLOCK
#RUN_BLOCK

#END_RUN_BLOCK
#POSTPROCESSING_BLOCK

#END_POSTPROCESSING_BLOCK
#PLOT_BLOCK

#END_PLOT_BLOCK
#CALCTRANSFORM_BLOCK

#END_CALCTRANSFORM_BLOCK
#SUBMITNEXT_BLOCK

#END_SUBMITNEXT_BLOCK
#BATCH_BLOCK

#END_BATCH_BLOCK
#CLEANUP_BLOCK

#END_CLEANUP_BLOCK''')



    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
           
####################################################################################################################


def __temp_executable(oneCalc,ID,from_step=0):
    '''
    Creates the _ID.py for a single calculation for a given step.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	  ID (str): ID string for the particular calculation and step
	  
    Keyword Argument
	  ext (str): an optional postfix to the ID's to all calculations in the set for a given step

    Returns:
         The one calculation for the new set, its ID, the old calculation, and its ID

    '''

    oneCalc['__chain_index__']+=1
    try:
	    pass
    except:
	    pass
    

    temp_calcs=copy.deepcopy(oneCalc)
    new_calcs=copy.deepcopy(oneCalc)



    logging.debug("ENTERING __temp_executable")
    try:

            f=ID
            d = temp_calcs

            subdir = temp_calcs['_AFLOWPI_FOLDER_']
            a = f+'.out'


            """generate a dummy .py file for use later when scf has run and we have info to generate hash"""
            finp = oneCalc['execFile']

	    if from_step==0:
		    extension= '_%.02d' % oneCalc['__chain_index__']
	    else:
		    extension= '_%.02d' % from_step
#	    extension= '_step_%.02d' % (int(ID.split('_')[2]))
	    
#	    if ext!='':
#		    extension+='_%s'%ext


            temp_hash = '%s%s'%(ID.split('_')[0],extension)
#	    if temp_hash==oneCalc['_AFLOWPI_PREFIX_'][1:]+'_step_01':
#		    temp_hash=temp_hash.split('_')[0]

		    

		    
            finp_temp = "_%s.py" % (temp_hash)
            finp = os.path.join(subdir,finp_temp)
            



	    
	    AFLOWpi.prep.__writeTemplate(finp)
            
            new_calcs['execFile']=finp
#	    print finp
#            oneCalc['nextCalc']=temp_hash

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    logging.debug("EXITING __temp_executable")
    
    return new_calcs,temp_hash,oneCalc,ID

####################################################################################################################

def __removeFromBlock(oneCalc,ID,block,removal):
    '''
    Searches for a regular expression pattern inside a specific block of a each calculation's _ID.py and 
    removes the command from the block.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	  ID (str): ID string for the particular calculation and step
	  block (str): string of the block in the _ID.py that the addition is to be added to
	                for that step of workflow
	  removal (str): a regular expression to find and delete a pattern of text in a single block
	                  of each calculation's _ID.py for a given step

    Returns:
          None

    '''

    logging.debug('entering __removeFromBlock')
    subdir = oneCalc['_AFLOWPI_FOLDER_']    
    fileName = '_'+ID+'.py'

    try:
        with open(os.path.join(subdir,fileName),'r') as inputfile:
            inputFileText = inputfile.read()

        isolatedString=re.split('#%s_BLOCK' % block,inputFileText)
        isolatedStringEnd = isolatedString[-1]
        isolatedStringBeginning = isolatedString[0]
        isolatedStringEmpty=re.split('#END_%s_BLOCK' % block,isolatedStringEnd)[-1]
        isolatedString=re.split('#END_%s_BLOCK' % block,isolatedStringEnd)[0]
    except:
        pass

    remainingText = re.sub(removal,'',isolatedString)
    newBlockText = '''#%s_BLOCK
%s
#END_%s_BLOCK''' % (block,remainingText,block)
    new_IDPY=isolatedStringBeginning+newBlockText+isolatedStringEmpty
    
    with open(os.path.join(subdir,fileName),'w') as outputfile:
        outputfile.write(new_IDPY)

    logging.debug('exiting __removeFromBlock')


def addToBlockWrapper(oneCalc,ID,block,addition):
	'''
	Wraps AFLOWpi.prep.__addToBlock for use inside _calcs_container methods

	Arguments:
              oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	      ID (str): ID string for the particular calculation and step
	      block (str): string of the block in the _ID.py that the addition is to be added to
	                    for that step of workflow
	      addition (str): a string containing code to be written to the specific block
	                       in the _ID.py for each calculation

        Returns:
              None
	
	'''

	AFLOWpi.prep.__addToBlock(oneCalc,ID,block,addition)

def __addToAll(calcs,block=None,addition=None):
    '''
    Adds text to a particular command block in the _ID.py for all calculations
    in the set.

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations 
      
    Keyword Arguments:
          block (str): the name of the command block that text will be added to
	  addition (str): the text to be added to the block
	  
    Returns:
          None

    '''

    if block==None or addition==None:
        return
    for ID,oneCalc in calcs.iteritems():
        __addToBlock(oneCalc,ID,block,addition)    



def __addToBlock(oneCalc,ID,block,addition):
    '''
    Writes the code to a single calculation's _ID.py

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	  ID (str): ID string for the particular calculation and step
	  block (str): string of the block in the _ID.py that the addition is to be added to
	                for that step of workflow
	  addition (str): a string containing code to be written to the specific block
	                   in the _ID.py for each calculation

    Returns:
          None

    '''
 
    logging.debug('entering __addToBlock')
    if '__DONOTBUILD__' in globals().keys():
        logging.debug('not writing to tree because build=False flag included in AFLOWpi.maketree')
    else:
        subdir = oneCalc['_AFLOWPI_FOLDER_']    
        fileName = '_'+ID+'.py'
        try:
            with open(os.path.join(subdir,fileName),'r') as inputfile:
                inputFileText = inputfile.read()
#                logging.debug(inputFileText)
                insert =  addition

                newText = r'%s\n#END_%s_BLOCK' % (insert,block)
                newinputfile  =  re.sub('#END_'+block+'_BLOCK',newText, inputFileText,re.MULTILINE)

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        try:
            with open(os.path.join(subdir,'_'+ID+'.py'),'w') as outfileObj:
                outfileObj.write(newinputfile)

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

    logging.debug('exiting __addToBlock')    
    ##########################################################################################################################

def __writeLoggingBlock(oneCalc,ID,logfile):
	'''
	###############
	###OBSOLETE###
	###############
	Writes the mechism to start the logging at runtime in the _ID.py

	Arguments:
     	      oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
 	      ID (str): ID string for the particular calculation and step
	      logfile (str): string of the path of the logfile within the directory tree for that set of calculations

	Returns:
	      None

	'''
	     
	__addToBlock(oneCalc,ID,'LOGGING','import logging')
	__addToBlock(oneCalc,ID,'LOGGING','logfile="%s"' % logfile)
	__addToBlock(oneCalc,ID,'LOGGING','''logging.basicConfig(filename=logfile,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',level=AFLOWpi.prep.__getLoglevel())''')



def __writeToScript(executable,oneCalc,ID,from_step=0):
	"""
	Generates calls on several functions to set up everything that is 
        needed for a new step in the workflow. AFLOWpi.prep.writeToScript loops
	over the calculation set and calls this function to do the work on
	each calculation in the set.

	Arguments:
     	      oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
 	      ID (str): ID string for the particular calculation and step
	      executable (str): <DEFUNCT OPTION: HERE FOR LEGACY SUPPORT> 
	      *args: <DEFUNCT OPTION: HERE FOR LEGACY SUPPORT>

        Keyword Arguments:
	      **kwargs: <DEFUNCT OPTION: HERE FOR LEGACY SUPPORT>

        Returns:
	      A single calculation dictionary for a new step and its ID.
	      
	"""

        logging.debug('entering __writeToScript')
        # try:
        #     if 'configFile' in kwargs.keys():
        #         if kwargs['configFile'] == None:
        #             configFile= __getConfigFile()
        #         else:
        #             configFile = kwargs['configFile']


        #     else:
        #         configFile=None
        # except Exception,e:
        #     AFLOWpi.run._fancy_error_log(e)
	# try:
	# 	new_job = kwargs['new_job']
	# 	del kwargs['new_job']
	# except:
	# 	pass


	
	# if index<1:
		
	# 	ID=ID.split('_')[0]
	# 	print ID
		
	# print ID
	# raise SystemExit

        oneCalcCopy=copy.deepcopy(oneCalc)
        logfile = os.path.abspath(os.path.join((os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])),'AFLOWpi','LOG.log'))
	
	calc_type=''
#	if 'extension' in kwargs.keys():
#		calc_type=kwargs['extension']
#		del kwargs['extension']


        new_oneCalc,new_ID,oneCalc,ID = __temp_executable(oneCalcCopy,ID,from_step=from_step)
	index = oneCalc['__chain_index__']
        # try:
        #     argsString=''
        #     for k,v in kwargs.iteritems():
        #         if type(v)==type(4):
        #             oneKwarg = ('%s=%s,' % (k,v))
        #         elif type(v)==type('a'):
        #             oneKwarg = ('%s="%s",' % (k,v))
        #         elif type(v)==type(0.1):
        #             oneKwarg = ('%s=%s,' % (k,v))
        #         elif type(v)==type(True):
        #             oneKwarg = ('%s=%s,' % (k,v))
        #         else:
        #             oneKwarg=''
        #         argsString+=oneKwarg
        
        # except Exception,e:
        #     AFLOWpi.run._fancy_error_log(e)

        try:
            nextIDName,nextCalcName = __getNextOneCalcVarName(oneCalc,ID)
#            __addToBlock(new_oneCalc,new_ID,'PREPROCESSING','oneCalc,ID = %s(oneCalc,ID,"%s",%s)' % (executable,new_ID,argsString))
#            __addToBlock(new_oneCalc,new_ID,'PREPROCESSING',tryLoadStr % (executable,new_ID,argsString))

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        try:

            


            prev_ID = '%s_%.02d' % (ID.split('_')[0],from_step-1)

#            prev_ID = '%s_step_%.02d' % (ID.split('_')[0],oneCalc['__chain_index__']-1)



            if os.path.exists(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+prev_ID+'.py')):
                pass

            else:
                pass
#                prev_ID=ID
#	    print prev_ID
#	    print ID
#	    print new_ID

	    oneCalc['__chain_index__']=from_step
#	    print oneCalc['__chain_index__']
            force_one_job=False
	    try:
		    if new_job:
			    force_one_job=False
		    else:
			    force_one_job=True
	    except:
		    new_job=False
	    __addToBlock(new_oneCalc,new_ID,'LOCK',"oneCalc,ID = AFLOWpi.prep.__lock_transform(oneCalc,ID)")            

	    
	    if new_ID!=oneCalc['_AFLOWPI_PREFIX_'][1:]:
                #add one to the index counter for the next step
                __addToBlock(oneCalc,prev_ID,'SUBMITNEXT',"globals()['__TEMP__INDEX__COUNTER__']+=1")                             
                __addToBlock(oneCalc,prev_ID,'SUBMITNEXT','''AFLOWpi.run.__submitJob("%s",oneCalc,__submitNodeName__,forceOneJob=%s)''' % (new_ID,force_one_job))
                __addToBlock(oneCalc,prev_ID,'SUBMITNEXT',"globals()['__TEMP__INDEX__COUNTER__']-=1")

            __fillTemplate(oneCalc,new_ID)
            '''at a global start timer for figuring out runtime of pw.x'''

            __addToBlock(new_oneCalc,new_ID,'LOADCALC','''oneCalc = AFLOWpi.prep.__loadOneCalc('%s','%s')''' % (oneCalc['_AFLOWPI_FOLDER_'],ID))

            tryLoadStr = '''
try:
    oneCalc = AFLOWpi.prep.__loadOneCalc('%s','%s')
except:
    pass

oneCalc['__chain_index__']=%d
''' % (oneCalc['_AFLOWPI_FOLDER_'],new_ID,from_step)
            __addToBlock(new_oneCalc,new_ID,'LOADCALC',tryLoadStr)
#            __addToBlock(oneCalc,ID,'CALCTRANSFORM','''oneCalc,ID = AFLOWpi.prep.__getEfermi(oneCalc,ID)''') 
            __addToBlock(oneCalc,ID,'CALCTRANSFORM','''AFLOWpi.prep.__saveOneCalc(oneCalc,ID)''') 
            __addToBlock(oneCalc,new_ID,'IMPORT','globals()["__TEMP__INDEX__COUNTER__"]=%s'%from_step)

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            

        try:
            with open(oneCalc['execFile'],'a') as prevExecFile:
                logging.debug('   ID in prep %s' % ID)
                configFile= __getConfigFile()

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        logging.debug(new_ID)
        logging.debug('exiting __writeToScript')            
        
        new_oneCalc['__calcVarName__']= nextCalcName
        new_oneCalc['__IDVarName__']= nextIDName
        new_oneCalc['__execCounterBkgrd__']=0
        new_oneCalc['__execCounter__']=0
#	try:
#		new_oneCalc['__chain_index__']+=1
#	except Exception,e:
#		new_oneCalc['__chain_index__']=1
#        new_oneCalc['__chain_index__']=globals()['__TEMP__INDEX__COUNTER__']
        
        return new_oneCalc,new_ID



#########################################################################################################################
#########################################################################################################################
## merge from multiple calcs
#########################################################################################################################
######################################################################################################################### 

	


def __loadCalcsFromOneCalc(oneCalc,ID):
    '''
    Loads the entire calc set from one of the calculations at runtime. Used when all jobs in 
    a set are needed to perform a task (i.e plotting data from all calculations in the set)

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
 	  ID (str): ID string for the particular calculation and step    

    Returns:
          The dictionary of dictionaries representing the calculation set that oneCalc 
	  is a part of

    '''
    
    PROJECT=oneCalc['PROJECT']
    SET=oneCalc['SET'] 
    
    CONFIG=os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'AFLOWpi','CONFIG.config')
    try:
	    #THIS NEEDS TO BE BETTER!!!!!
        logname='step'+re.findall(r'_\d+',ID)[-1]
    except:
        logname='step_01'

    calcs = AFLOWpi.prep.loadlogs(PROJECT,SET,logname,config=CONFIG)

    return calcs

def __checkSuccessCompletion(oneCalc,ID,faultTolerant=True):
    '''
    Allows one of the calculations to check the status of the other calculations
    and returns True if all calculations in the set have completed.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
 	  ID (str): ID string for the particular calculation and step    
    

    Keyword Arguments:
          faultTolerant (bool): a flag to choose if we return True if some of the
                                 calculations ran but did not complete successfully
	                     

    Returns:
          True if all calculations in the set are complete and False if not.

    '''

    calcs = __loadCalcsFromOneCalc(oneCalc,ID)



    for ID,oneCalc in calcs.iteritems():
        logging.debug('status of %s: %s is: %s' % (oneCalc['_AFLOWPI_INDEX_'],ID,oneCalc['__status__']['Complete']))
        if oneCalc['__status__']['Complete']==False:
            if faultTolerant==True:
                '''if a calc failed for some reason ignore the fact that it didn't complete'''
                if oneCalc['__status__']['Error']!='None':
                    pass
                else:
                    return False

            else:
                '''if fault tolerance is false then return false if any of the calcs did not complete'''
                return False


    '''if all calcs are completed for this step return True'''
    return True


	##########################################################################################
	###NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED ###
	##########################################################################################
# def __saveMasterList(calcs,__master__key__):

#     baseDir = os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])
    

#     AFLOWpiDir=os.path.join(baseDir,'AFLOWpi')

# def __loadMasterList(oneCalc,ID,__master__key__):
#     baseDir = os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])
#     AFLOWpiDir=os.path.join(baseDir,'AFLOWpi')
#     masterEXEList = os.path.join(AFLOWpiDir,'%s.exeList' %__master__key__)


def runAfterAllDone(calcs,command,faultTolerant=True):
    '''   
    Adds a command to the BATCH command block at the end of each calculation's _ID.py 
    for all calculations in the set. Used to execute a command over all calculations
    in particular step have completed.

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations 

    Keyword Arguments:
	  command (str): the text to be added to the BATCH block
          faultTolerant (bool): a flag to choose if we return True if some of the
                                 calculations ran but did not complete successful          
    
    Returns:
          None

    '''

    completionBool='''if AFLOWpi.prep.__checkSuccessCompletion(oneCalc,ID,faultTolerant=%s):''' % faultTolerant
    addition = """
import copy
AFLOWpi.prep.__saveOneCalc(oneCalc,ID)


%s
         try:
            calcs = AFLOWpi.prep.__loadCalcsFromOneCalc(oneCalc,ID)
         except Exception,e:
            AFLOWpi.run._fancy_error_log(e)



         %s
else: 
         pass

""" % (completionBool,command)
    __addToAll(calcs,block='BATCH',addition=addition)
    


#########################################################################################################################
#########################################################################################################################
## hashing
#########################################################################################################################
######################################################################################################################### 

import hashlib
def __hash64String(hashString):
    '''
    Takes a string and returns a 16 character long CRC64 hash checksum of it.

    Arguments:
          hashString (str): the string you want a hash checksum of


    Returns:
          hashed_str (str): a 16 character long CRC64 hash of the inputted string

    '''

    hashString = __parseRef(hashString)
    return __crc64digest(hashString)


def __crc64digest(aString):
    """
    Generates the CRC64 hash checksum of the inputted string.

    From W. H. Press, S. A. Teukolsky, W. T. Vetterling, and 
    B. P. Flannery, "Numerical recipes in C

    Arguments:
          aString (str): a string

    Returns:
          hashed_str (str): a 16 character long CRC64 hash of the inputted string
    """

    CRCTableh = [0] * 256
    CRCTablel = [0] * 256
    def _inittables(CRCTableh, CRCTablel, POLY64REVh, BIT_TOGGLE):
            for i in xrange(256):
                    partl = i
                    parth = 0L
                    for j in xrange(8):
                            rflag = partl & 1L
                            partl >>= 1L
                            if parth & 1:
                                    partl ^= BIT_TOGLE
                            parth >>= 1L
                            if rflag:
                                    parth ^= POLY64REVh
                    CRCTableh[i] = parth
                    CRCTablel[i] = partl
    # first 32 bits of generator polynomial for CRC64 (the 32 lower bits are
    # assumed to be zero) and bit-toggle mask used in _inittables
    POLY64REVh = 0xd8000000L
    BIT_TOGGLE = 1L << 31L
    # run the function to prepare the tables
    _inittables(CRCTableh, CRCTablel, POLY64REVh, BIT_TOGGLE)
    # remove all names we don't need any more, including the function
    del _inittables, POLY64REVh, BIT_TOGGLE
    # this module exposes the following two functions: crc64, crc64digest
    def crc64(bytes, (crch, crcl)=(0,0)):
            for byte in range(len(bytes)):
                    shr = (crch & 0xFF) << 24
                    temp1h = crch >> 8L
                    temp1l = (crcl >> 8L) | shr
                    tableindex = (crcl ^ ord(bytes[byte])) & 0xFF 
                    crch = temp1h ^ CRCTableh[tableindex]
                    crcl = temp1l ^ CRCTablel[tableindex]
            return crch, crcl


    checksum="%08X%08X" % (crc64(aString))
    return checksum.lower()


def cleanCalcs(calcs,runlocal=False):
    '''
    Wrapper function for AFLOWpi.prep.__cleanCalcs

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations     

    Keyword Arguments:
	  runlocal (bool): a flag to choose whether or not to run the wrapped function now
	                    or write it to the _ID.py to run during the workflow   

    Returns:
          None

    '''

    try:
        for ID,oneCalc in calcs.iteritems():
            if runlocal:
                AFLOWpi.retr.__cleanCalcs(ID,oneCalc)
            else:
                AFLOWpi.prep.__addToBlock(oneCalc,ID,'CLEANUP','AFLOWpi.prep.__cleanCalcs(ID,oneCalc)')
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

def __cleanCalcs(ID,oneCalc):
    '''
    Removes all files from the directory tree with a "_" prefix.
    
    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step
    
    Returns:
          None

    '''

    fileList = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/_*')
    for item in fileList:
        try:
            os.remove(item)
        except OSError:
            shutil.rmtree(item)


#########################################################################################################################
#########################################################################################################################
## calc transforms
#########################################################################################################################
######################################################################################################################### 


def bands(calcs,dk=None,nk=None):
	'''
	Wrapper function to write the function AFLOWpi.prep.__oneBands to the _ID.py

	Arguments:
	      calcs (dict): a dictionary of dicionaries representing the set of calculations 

	Keyword Arguments:
	      dk (float): distance between points for Electronic Band Structure calculation
	      nk (int): approximate number of k points to be calculated along the path
		  	      
	Returns:
	      The identical "calcs" input variable

	'''

	loadModString = 'AFLOWpi.prep.__oneBands(oneCalc,ID,dk=%s,nk=%s)'%(dk,nk)
	__addToAll(calcs,block='PREPROCESSING',addition=loadModString)
	return calcs



@newstepWrapper(__check_lock)
def __oneBands(oneCalc,ID,dk=None,nk=None,configFile=None):
	"""
	Add the ncsf with the electronic band structure path k points input to each 
	subdir and update the master dictionary.

	Arguments:
     	      oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
 	      ID (str): ID string for the particular calculation and step

	Keyword Arguments:
	      dk (float): the density in the Brillouin zone of the k point sampling along the 
                           entirety of the path between high symmetry points.
	      nk (int): the approximate number of sampling points in the Brillouin Zone along
 	                 the entirety of the path between high symmetry points. Points are chosen 
		         so that they are equidistant along the entirety of the path. The actual 
		         number of points will be slightly different than the inputted value of nk.
		         nk!=None will override any value for dk.


	Returns:
	      A calculation that has been transformed from scf to nscf with k point sampling
	      along the high symmetry path in the Brillouin Zone

	"""

        logging.debug('entering __oneBands')
        try:
                inputfile = oneCalc['_AFLOWPI_INPUT_']                
                '''we need nscf and not bands because we need output for HOMO'''
                subdir = oneCalc['_AFLOWPI_FOLDER_']

                if not os.path.exists(os.path.join(subdir,'%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:])):
                    logging.error('%s does not exist. Can not make bands calculations. Exiting' % os.path.join(subdir,'%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:]))
                    raise SystemExit


        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


        '''checks to see if "crystal_b" has already been substituted for "automatic"'''
        try:
            prevID_str=oneCalc['prev'][-1]
            for prev_calc in  oneCalc['prev']:
                try:
                    outfile = open(os.path.join(subdir,'%s.out' % prev_calc),'r').read()
                    match1 = re.findall(r'Kohn-Sham states\s*=\s*([0-9]*)',outfile)[-1]
                    nbnd = int(2.0*int(match1))
                    prevID_str=prev_calc
                    break
                except:
			pass

            bandsKRegex = AFLOWpi.qe.regex.k_points('','content','regex')
            '''sets dk to generic value if no option for dk or nk is used'''
            if nk==None and dk==None:
                dk=0.1
                '''if dk isn't set but nk is..set dk small and average out later with nk'''
            elif dk==None and nk!=None:
                dk=0.001


            path = AFLOWpi.retr.__getPath(dk,oneCalc,ID=prevID_str)


            if nk!=None:
                '''scale the k points in the path list so they're as close to nk as possible'''
                splitPath =  [[x.split()[0],int(x.split()[1])] for x in  path.split('\n')[1:] if len(x)!=0]
                total =  [int(x.split()[1]) for x in  path.split('\n')[1:] if len(x)!=0]

                for entries in range(len(total)):
                    if total[entries]==0:
                        total[entries]+=1

                total=sum(total)
                scaleFactor = float(nk)/total



		path = AFLOWpi.retr.__getPath(dk/scaleFactor,oneCalc,ID=prevID_str)


            inputSplit=AFLOWpi.retr.__splitInput(inputfile)


            inputSplit['&control']['calculation']="'nscf'"
            inputSplit['&system']['noinv']=".true."
            inputSplit['&system']['nosym']=".true."
	    try:
		    inputSplit['&system']['nbnd']=nbnd
	    except:
		    pass

            inputSplit['K_POINTS']['__modifier__']='{crystal_b}'
            inputSplit['K_POINTS']['__content__']=path
#            inputSplit['&control']['wf_collect']='.TRUE.'
	    
            inputfile = AFLOWpi.retr.__joinInput(inputSplit)

            calc_label = ID

            with open(os.path.join(subdir,ID+'.in'),'w') as new_inputfile:
                new_inputfile.write(inputfile)


            inputSplit=AFLOWpi.retr.__splitInput(inputfile)
            inputfile = AFLOWpi.retr.__joinInput(inputSplit)

###############################################################################################################


        #     a1 = '_'+temp_ID+'.py'

	    oneCalc['_AFLOWPI_INPUT_'] = inputfile
        #     d['execFile']= os.path.join(subdir,'_%s.py' % temp_ID)

        #     try:
        #         __fillTemplate(d,temp_ID)
        #     except Exception,e:
        #         AFLOWpi.run._fancy_error_log(e)

        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            logging.error("%s not in %s" % (ID,oneCalc['_AFLOWPI_FOLDER_']))
            
        # logging.debug('exiting __oneBands')

        # '''set prev to current ID so that we can check if we have already transformed'''
        # d['prev'].append(temp_ID)
        
        # '''update status of the calc after it's transformed'''
        # d['__status__']=collections.OrderedDict({"Start":False,'Complete':False,"Restart":0,"Error":'None'})

        # """save it so if it's loaded by a scipt it'll have the correct vals"""
        # __saveOneCalc(d,temp_ID)

	return oneCalc,ID
####################################################################################################################
import math

def bandsAflow(dk, LAT):
	"""
	Query aflow for band structure path and generate the path for band structure calculation

	Arguments:
	      dk (float):  distance between k points along path in Brillouin Zone
	      LAT (int): bravais lattice number from Quantum Espresso convention

	Keyword Arguments:
	      None

	Returns:
	     info (str): some information
	     nks (str): number of k points in path
	     stringk (str): kpoint path string
	"""

	#logging.info('Running aqe')
	COMMAND = ['./aqe', '--bzd', str(LAT)]
	proc = subprocess.Popen(COMMAND, stdout=subprocess.PIPE)
	fromaqe = proc.communicate()[0]	

# we use aqe --qe --bzd LAT
# e.g. ./aqe --bzd RHL1
#RHL1 (rhombohedral alpha<90) G-L-B1 B-Z-G-X Q-F-P1-Z L-P
#20   ! 20 grids 
#Line-mode
#reciprocal
#   0.000   0.000   0.000    ! \Gamma
#   0.500   0.000   0.000    ! L
#
#   0.500   0.000   0.000    ! L
#   0.500  0.5  -0.5  ! B_1
#
#   0.5  0.500  0.5  ! B
#   0.500   0.500   0.500    ! Z
#
#   0.500   0.500   0.500    ! Z
#   0.000   0.000   0.000    ! \Gamma
#
#   0.000   0.000   0.000    ! \Gamma
#   0.5  0.000  -0.5  ! X
#
#   0.5  0.5  0.000   ! Q
#   0.500   0.500   0.000    ! F
#
#   0.500   0.500   0.000    ! F
#   0.5  0.5  0.5  ! P_1
#
#   0.5  0.5  0.5  ! P_1
#   0.500   0.500   0.500    ! Z
#
#   0.500   0.000   0.000    ! L
#   0.5  0.5  0.5  ! P
	
	a = re.split('\n',fromaqe)
	kkk = {}
	info = re.split('\s',a[0])
	for l in a[4:]:
		if l != '':
			b = re.split('\s+',l)
			label = re.sub(r'\\|_|amma','', b[-1])
			kk = []
			for k in b:
				try:
					kk.append(float(k))
				except:
					pass
			kkk.update({label: tuple(kk)})
	
	path = []
	stringk = ''
	for l in info:
		if re.search('\w-\w',l):
			p = re.split('-',l)
			path.append(p)

	nk = []
	for piece in path:
		ks = kkk[piece[0]]
		for k in piece:
			ke = kkk[k] 
			km = 0.0
			for i in range(3): km  += (ke[i]-ks[i])**2
			nk.append(int(math.sqrt(km)/dk))
			ks = ke

	
	jc = 1
	for piece in path:
		for k in piece:
			t1 = '% 05.3f % 05.3f % 05.3f' % kkk[k]
			try:
				t2 = ' %3d' % nk[jc]
			except:
				t2 = ' %3d' % 0
			t3 = ' !'+k+'\n'
			stringk += t1+t2+t3
			jc += 1
	nks = jc - 1	
	#logging.debug('Exiting bandsAflow')
		
	return info, nks, stringk


import filecmp
#import AFLOWpi
import os
import shutil




def maketree(calcs, pseudodir=None,build=True,workdir=None):
	"""
	Make the directoy tree and place in the input file there
	
	Arguments:
	      calcs (dict): - Dictionary of dictionaries of calculations

	Keyword Arguments:
	      pseudodir (str): path of pseudopotential files directory 
              workdir (str): a string of the workdir path that be used to override what is in the 
	         	   config file used when initating the AFLOWpi session
	      build (bool): <DEFUNCT OPTION. NEEDS REMOVAL>


        Returns:
	      None

	"""

        logging.debug('Entering makeAFLOWpitree')
        megahashStr=''
        for ID in calcs.keys():
            megahashStr+=ID
            __MASTER__KEY__ = AFLOWpi.prep.__crc64digest(megahashStr)


        
        try:
            if pseudodir==None:
                pseudodir = __ConfigSectionMap('prep','pseudodir')
                if os.path.isabs(pseudodir) == False:
                    configFileLocation = __getConfigFile()
                    configFileLocation = os.path.dirname(configFileLocation)
                    pseudodir =  os.path.normpath(os.path.join(configFileLocation, pseudodir))

            else:
                if os.path.isabs(pseudodir) == False:
                    pseudodir =  os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), pseudodir))

            configFileName = __getConfigFile()

            topdir=''
            for ID,oneCalc in calcs.iteritems():
                topdir = os.path.abspath(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])))
                break
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

        printFlag=1
        if len(calcs)>10:
                printFlag=0
                print "Writing input files to directory tree"



        counter=0
        for k,v in calcs.iteritems():
            try:
                    # create the subdir under topdir
                    counter += 1

                    ddir = v['_AFLOWPI_FOLDER_']

                    logfile = os.path.abspath(os.path.join(os.path.dirname(ddir),'AFLOWpi','LOG.log')) 

                    if not os.path.exists(ddir):
                            os.mkdir(ddir)


                    """Get name of logfile for this calc"""
            except Exception,e:
                AFLOWpi.run._fancy_error_log(e) 

        if build==False:
            global __DONOTBUILD__
            __DONOTBUILD__= True
            
        else:

            for k,v in calcs.iteritems():
                try:

                        ddir=v["_AFLOWPI_FOLDER_"]

                        """place the input file in the dir"""
                        finp = os.path.join(ddir,k)+'.in'

                        f = open(finp,'w')
                        f.write(v['_AFLOWPI_INPUT_'])
                        f.close()

                        """generate initial SCF executable in folder"""
                        calcs[k]['_AFLOWPI_CONFIG_'] = configFileName

                        calcs[k]['_AFLOWPI_FOLDER_']=ddir

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e) 

                try:
                    finp = os.path.join(ddir,'_'+k)+'.py'
                    calcs[k]['execFile']=finp                        
                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e) 

                try:

                    if '__DONOTBUILD__' in globals().keys():
                        logging.debug('not writing to tree because build=False flag included in AFLOWpi.maketree')
                    else:

			    __writeTemplate(finp)
			    __fillTemplate(v,k)


                    configFile= __getConfigFile()

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e) 

                try:
                    clusterType = AFLOWpi.prep.__ConfigSectionMap("cluster",'type').upper()
                    if clusterType == 'PBS':
                        qsubtext = AFLOWpi.run.__qsubGen(v,k)

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)


                # copy the pseudo file
                if pseudodir != None:
                        for l,m in v.iteritems():
                                if re.search(r'_AFLOWPI_[A-Z]\d*PSEUDO_',l):
                                    try:
                                                a = os.path.join(pseudodir,m)
                                                if __ConfigSectionMap('prep','copypseudos').lower() == 'false':
                                                    pseudodir = __ConfigSectionMap('prep','pseudodir')
                                                    if os.path.isabs(pseudodir) == False:
                                                        configFileLocation = __getConfigFile()
                                                        configFileLocation = os.path.dirname(configFileLocation)
                                                        pseudodir =  os.path.join(configFileLocation, pseudodir)


                                                    if os.path.exists(os.path.join(pseudodir,m)):
                                                        try:
                                                            os.symlink(os.path.join(pseudodir,m),os.path.join(v['_AFLOWPI_FOLDER_'],m))
                                                        except OSError,e:
                                                            try:
                                                                os.unlink(os.path.join(v['_AFLOWPI_FOLDER_'],m))
                                                                os.symlink(os.path.join(pseudodir,m),os.path.join(v['_AFLOWPI_FOLDER_'],m))
                                                            except:
                                                                print m
                                                        except Exception,e:
                                                            print e
                                                    else:
                                                        logging.error('Count not find %s check if your pseudodir in config is properly set' % os.path.join(pseudodir,m))
                                                    
                                                else:
                                                    if printFlag:
							    pass
#                                                            print 'Updating ', l
#                                                            print 'Copying '+a+' to '+ddir
                                                    if not os.path.exists(os.path.join(ddir,m)):
                                                        shutil.copy(a, ddir) 
                                    except IOError as e:
                                        print e
                                        logging.critical(e)
                                    except AttributeError as e:
                                        print e
                                        logging.critical(e)	

                                    
        '''saving the first set of calcs so they can be loaded by the script when the first calcs start'''
        for ID,oneCalc in calcs.iteritems():
            __saveOneCalc(oneCalc,ID)

            if __findInBlock(oneCalc,ID,'ONECALC','''oneCalc = AFLOWpi.prep.__loadOneCalc('%s','%s')''' % (oneCalc['_AFLOWPI_FOLDER_'],ID))==False:
                __addToBlock(oneCalc,ID,'ONECALC','''oneCalc = AFLOWpi.prep.__loadOneCalc('%s','%s')''' % (oneCalc['_AFLOWPI_FOLDER_'],ID))
                

                


#            if __findInBlock(oneCalc,ID,"ID","""__master__key__ = '%s'""" % __MASTER__KEY__)==False:
#                __addToBlock(oneCalc,ID,"ID","""__master__key__ = '%s'""" % __MASTER__KEY__)


#        try:
#            index= globals()['__TEMP__INDEX__COUNTER__']
#        except:
	index=1
	

        updatelogs(calcs,'step_%02d'%index)
        logging.debug('Exiting makeAFLOWpitree')


def totree(tobecopied, calcs,rename=None,symlink=False):
	"""
	Populate all the subdirectories for the calculation with the file in input

	Arguments:
	      tobecopied (str): filepath to be copied to the AFLOWpi directory tree
 	      calcs (dict): Dictionary of dictionaries of calculations
 
        Keyword Arguments:
	      rename (bool): option to rename the file/directory being moves into the AFLOWpi
	                  directory tree
	      symlink (bool): whether to copy the data to the AFLOWpi directory tree or 
	                   to use symbolic links

        Returns:
	      None
		
	"""

	logging.debug("Entering totree")

        try:
            topdir = os.path.abspath(os.path.dirname(calcs[random.choice(calcs.keys())]['_AFLOWPI_FOLDER_']))
            
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            
	printBool=1
	if len(calcs)>10:
		printBool=0
		print 'Copying '+tobecopied+' to directory tree'

        try:
            actualFile = os.path.split(tobecopied)[-1]
        except:
            pass
	for l,v in calcs.iteritems():

		a = v.get('_AFLOWPI_FOLDER_')
                if rename != None:
                    a = os.path.join(a,rename)
                else:
                    a = os.path.join(a,actualFile)
                if symlink==True:
                    try:
                        os.symlink(tobecopied,a)
                    except OSError,e:
                        pass
                else:
                    if printBool:
#                            print 'Updating ', l
#                           print 'Copying '+tobecopied+' to '+a
			    pass

                    try:

                        if os.path.isfile(tobecopied): 
                            try:
                                if os.path.exists(os.path.join(a,actualFile)):
                                    try:
					    if not filecmp.cmp(tobecopied,os.path.join(a,actualFile)):
						    shutil.copy(tobecopied, a) 
					
					    else:
						    if printBool:
							    logging.debug('Identical copy of %s Exists in %s, not copying to this subdirectory of tree' % (actualFile,a))
				    except:
					    pass
                                else:
					try:
						shutil.copy(tobecopied, a) 
					except:
						pass
                            except Exception,e:

				    AFLOWpi.run._fancy_error_log(e)

                    except IOError as e:
                            print tobecopied+" does not exist in AFLOWpi Directory"
                            logging.error(tobecopied+" does not exist in AFLOWpi Directory")
                            logging.error(e)
                    except Exception,e:
                        AFLOWpi.run._fancy_error_log(e)


	logging.debug("Exiting totree")



####################################################################################################################

####################################################################################################################





#########################################################################################################################
#########################################################################################################################
## CALCLOGS
#########################################################################################################################
#########################################################################################################################


def updatelogs(calcs,logname,runlocal=False):

    if logname=='LOG' or logname=='AFLOWKEYS' or logname=='summary':
        logging.warning('calc log name cannot be the same as master log file i.e. LOG,AFLOWKEYS, or summary')
        print 'calc log name cannot be the same as master log file i.e LOG,AFLOWKEYS, or summary'
        return



    '''find calclogdir and if it doesn't exist create it (it should exist)'''
    for ID,oneCalc in calcs.iteritems():
        baseDir=os.path.abspath(os.path.join((os.path.dirname(oneCalc['_AFLOWPI_FOLDER_'])),'AFLOWpi'))
        calcLogDir=os.path.join(baseDir,'calclogs')
        if not os.path.exists(calcLogDir):
            os.mkdir(calcLogDir)
        break
    '''write the location of each of the _ID.oneCalc files for all the calcs'''
    with open(os.path.join(calcLogDir,logname+'.log'),'w') as testLogFile:
        for ID,oneCalc in calcs.iteritems():
            testLogFile.write(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.oneCalc' % ID)+'\n')

        


def loadlogs(PROJECT='',SET='',logname='',config=None):

         print 'Reading from Config: %s' % config

         configFile=config
         if os.path.exists(config)==False:
             logging.info('configuration file: %s does not exist.' % config)
             print 'configuration file: %s does not exist.' % config
             raise SystemExit


         '''
         put configFile global into the namespace as it would be if we did AFLOWpi.prep.init
         '''
         AFLOWpi.prep.__forceGlobalConfigFile(configFile)        

         workdir = __ConfigSectionMap('prep','workdir')



         if os.path.isabs(workdir) == False:
             configFileLocation = __getConfigFile()
             configFileLocation = os.path.dirname(configFileLocation)
             workdir =  os.path.join(configFileLocation, workdir)

         """
         FIX THSI         FIX THIS         FIX THIS         FIX THIS         FIX THIS
         """
         if os.path.isabs(config) == False:
             execDIR=os.path.dirname(os.path.abspath(__main__.__file__))
             config = os.path.abspath(os.path.join(execDIR,config))
         else:
             pass
         
         AFLOWpi.prep.__forceGlobalConfigFile(config)
         """
         FIX THSI         FIX THIS         FIX THIS         FIX THIS         FIX THIS
         """
         filename = os.path.join(workdir,PROJECT,SET,'AFLOWpi','calclogs',logname+'.log')

         filename1= os.path.join(workdir,PROJECT,SET,'AFLOWpi','calclogs',logname+'.log')

         if os.path.exists(filename)==False:
             """in case we have split files with .dat,.bak,and .dir"""

             if os.path.exists(os.path.join(workdir,PROJECT,SET,'AFLOWpi','calclogs',logname+'.dat'))==False:

                 if os.path.exists(os.path.join(workdir,PROJECT,SET,'AFLOWpi','calclogs',logname+'.log'))==False:
                     
                     logging.info('saved calc log: %s does not exist. Check the workdir parameter in your config file as well as project and set of these calcs.' % filename)
                     print 'saved calc log: %s does not exist. Check the workdir parameter in your config file  as well as project and set of these calcs.' % filename
                     raise SystemExit
             else:
                 filename = os.path.join(workdir,PROJECT,SET,'AFLOWpi','calclogs',logname)
                 
         

         '''
         pull logs from the file
         '''
         with open(filename1,'r') as logLocFile:
             logLocString = logLocFile.read()
         logFiles=[x for x in logLocString.split('\n') if len(x)!=0]

         returnCalcs=collections.OrderedDict()
         for log in logFiles:
             ID = '.'.join(os.path.basename(log).split('.')[:-1])[1:]
             folder = os.path.dirname(log)
             '''
             try to load the log from each folder but if it doesn't exist
             or none of the ID.oneCalc files exist return a blank string
             '''
             try:
		     with open(log,'rb') as oneCalcPickleFile:
			     oneCalc = cPickle.load(oneCalcPickleFile)    
		     #oneCalc=__loadOneCalc(folder,ID)

		     returnCalcs[ID]=oneCalc
             except Exception,e:
                 pass


         try:



             '''
             set the global counter number to what it was when the calcs saved so we don't
             overwrite the previous steps when the counter starts from 1 again
             '''
#             for ID,oneCalc in returnCalcs.iteritems():
#                 globals()['__TEMP__INDEX__COUNTER__'] = oneCalc['__TEMP__INDEX__COUNTER__']
                 
#                 break
             '''
             when the calc logs are loaded it's assumed that the previous steps have been 
             completed so it assumes we're starting from this step and sets the step loaded's 
             ID in the _<prefix>.Walltime to stepsasjobs='true' so that the timing works
             correctly for the restart in the case of stepsasjobs='false' in the config
             '''

#             for ID,oneCalc in returnCalcs.iteritems():
#                 walltimeLog = __readWalltimeLog('%s.Walltime' % oneCalc['_AFLOWPI_PREFIX_'])
#                 try:
#                     output_wall = shelve.open('%s.Walltime' % oneCalc['_AFLOWPI_PREFIX_'],writeback=True)   
#                     output_wall[ID]['stepsasjobs']='true'
#                     output.close()
#                 except Exception,e:
#                     AFLOWpi.run._fancy_error_log(e)


             returnCalcs = collections.OrderedDict(sorted(returnCalcs.iteritems(), key=lambda x: x[1]['_AFLOWPI_INDEX_']))
             return returnCalcs
         except:
             return returnCalcs


def __updatecalclogs(calcs,inc=True):
	"""
	Save the log for aflowkeys

	Arguments:
	       aflowkeys (dict): keys from aflow
	       calcs (dict): Dictionary of dictionaries of calculations
	 Keyword Arguments:
    	       inc (bool): choose to increment the logs or not (True,False)
	"""

        if not os.path.exists(os.path.join(AFLOWpidir,fc)):
            inc=False
        if not os.path.exists(os.path.join(AFLOWpidir,fa)):
            inc=False


	baseDir = os.path.dirname(calcs[random.choice(calcs.keys())]['_AFLOWPI_FOLDER_'])
 	AFLOWpidir = os.path.abspath(os.path.join(os.path.dirname(calcs[random.choice(calcs.keys())]['_AFLOWPI_FOLDER_']),'AFLOWpi'))
	
	fc = 'CALCSMASTER.log'
	
	if inc:
		inccalc = []


		for l in os.listdir(AFLOWpidir):
			if re.search(fc+'\d+',l):
				inccalc.append(int(l[15:]))



		if not len(inccalc):
			mc = 1
		else:
			mc = max(inccalc)+1

		fou = '%s%02d' % (fc, mc)
		A = os.path.join(AFLOWpidir, fc)
		B = os.path.join(AFLOWpidir, fou)
		os.rename(A,B)
		output = open(A,'w')
		cPickle.dump(calcs,output)
		output.close()
	else:
		output = open(os.path.join(AFLOWpidir,fc),'w')
		cPickle.dump(calcs,output)
		output.close()








import socket


def __saveOneCalc(oneCalc,ID):
#	print AFLOWpi.run.__get_index_from_pp_step(oneCalc,ID)
#	if ID==AFLOWpi.run.__get_index_from_pp_step(oneCalc,ID):

		with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.oneCalc' % ID),'wb') as oneCalcPickleFile:
			oldOneCalc = cPickle.dump(oneCalc,oneCalcPickleFile)    

def __loadOneCalc(folder,ID):
#	try:
#		ID=AFLOWpi.run.__get_index_from_pp_step(oneCalc,ID)
#	except:
#		pass
	with open(os.path.join(folder,'_%s.oneCalc' % ID),'rb') as oneCalcPickleFile:
		oldOneCalc = cPickle.load(oneCalcPickleFile)    
        try:
            with open(os.path.join(folder,'%s.in' % ID),'r') as inputFileObj:
                inputFileString = inputFileObj.read()
            
            oldOneCalc['_AFLOWPI_INPUT_']=inputFileString
            AFLOWpi.prep.__saveOneCalc(oldOneCalc,ID)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)

	return oldOneCalc
	



def __updatelogs(ID,oneCalc,filename):
    logging.debug('entering __updatelogs')
    """
    get index of this calc we're saving in the chain so we can
    restart with it if we do loadlogs in another script.
    """
    index=1
    try:
	    index = oneCalc['__chain_index__']
    except:
	    index =int(re.findall('step_(\d*)_*',ID)[-1])


#    oneCalc['__TEMP__INDEX__COUNTER__']=index
    """
    save this calc just in case it's the last in the current chain and we want to restart from it
    """
    __saveOneCalc(oneCalc,ID)



#########################################################################################################################
#########################################################################################################################
## LOCAL SCRATCH
#########################################################################################################################
#########################################################################################################################

def __scp_wrapper(from_path,to_path,node=''):


    try:

         '''try rsync first and fallback to scp'''
	 command="ssh -o StrictHostKeyChecking=no %s 'rsync -ru  %s %s/ > /dev/null' > /dev/null"%(node,from_path,to_path)
	 


	 if node=='':
		 command='cp -r  %s %s/'%(from_path,to_path)    
	 
	 try:
		 start=time.time()
		 os.system(command)
		 end=time.time()
	 except:
		 logging.debug('rsync failed. trying scp.')
		 command="ssh -o StrictHostKeyChecking=no %s 'cp -r %s %s/ > /dev/null' > /dev/null"%(node,from_path,to_path) 
		 os.system(command)



    except Exception,e:
	    pass


    diff=end-start
    size=0.0
    rate=0.0
    try:
	    
	    size = float(os.path.getsize(from_path))/(1024.0**2) 
	    if size<0.000001:
		    try:
			    size = float(os.path.getsize(to_path))/(1024.0**2) 
		    except:
			    pass

    except:
	 size=0.0
    try:
	 rate=size/diff
    except:
	 rate=0.0        

    logging.debug('\nfrom: %s\nto  : %s\nnode: %s\nsize: %s MB\ntime: %s sec\nrate: %s MB/s\n' %(from_path,to_path,node,size,diff,rate))



import copy
def __from_local_scratch(oneCalc,ID,ext_list=['save','wfc','hub','mix','occup','update','bfgs','restart','restart_k','restart_scf','ewfcp','ewfc','ewfcm']):
#def __from_local_scratch(oneCalc,ID,ext_list=['wfc','hub','mix','update','bfgs']):

    if AFLOWpi.prep.__ConfigSectionMap('cluster','localscratch').lower()=='true':      
        orig_list=copy.deepcopy(ext_list)
        try:
            try:
		    temp_dir= os.environ['TMPDIR']
	    except:
		    return
            with open(os.environ['PBS_NODEFILE'],'r') as nf:
                nf_string = nf.read()
            nf_list = list(set([x for x in nf_string.split('\n') if len(x.strip())!=0]))
            logging.debug(nf_list)

            work_dir=oneCalc['_AFLOWPI_FOLDER_']
#            ext_list=['wfc','hub','igk','restart','save','mix']
    
            prefix = AFLOWpi.retr.__prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])            

            ext_list=['%s.%s' % (prefix,x) for x in ext_list]


        except Exception,e:
            logging.debug(e)

        try:
            proc_by_node = [x for x in nf_string.split('\n') if len(x.strip())!=0]
            num_procs=len(proc_by_node)



            for ext in range(len(ext_list)):
                logging.info('copying %s files from:\nfrom : %s\nto   : %s' % (orig_list[ext],temp_dir,work_dir))
                proc_pool = multiprocessing.Pool(processes=num_procs)
                for proc in range(len(proc_by_node)):

                    try:
                        try:
				if orig_list[ext] in ['save','occup','update','bgfs']:
					file_name=os.path.join(temp_dir,'%s'%(ext_list[ext]))
				elif orig_list[ext] in ['restart','restart_k','igk','restart_scf']:
					if proc==0:
						file_name=os.path.join(temp_dir,'%s'%ext_list[ext])
					else:
						file_name=os.path.join(temp_dir,'%s%s'%(ext_list[ext],proc+1))

				else:
					file_name=os.path.join(temp_dir,'%s%s'%(ext_list[ext],proc+1))
                        except Exception,e:
                            continue


                        #add this transfer proc to the pool
			"""assuming the first file is in the folder we're executing this from"""
			if not os.path.exists(file_name) and proc==0:
				'''if we don't find the first skip to the next type'''
				break

                        res = proc_pool.apply_async(AFLOWpi.prep.__scp_wrapper, (file_name,work_dir,proc_by_node[proc])) 
			if orig_list[ext] in ['save','occup','update','bgfs',]:
				break

                    except Exception,e:
                        logging.debug(e)
                        pass

                #closes pool for this specific ext type and joins to wait for all procs to finish
                proc_pool.close()
                proc_pool.join()


        except Exception,e:
            logging.debug(e)
            pass




def __to_local_scratch(oneCalc,ID,ext_list=['save','wfc','hub','mix','occup','update','bgfs','restart','restart_k','restart_scf','ewfcp','ewfc','ewfcm']):
#def __to_local_scratch(oneCalc,ID,ext_list=['wfc','hub','mix','update','bgfs']):
    localscratchOpt=AFLOWpi.prep.__ConfigSectionMap('cluster','localscratch')
    localscratchOpt=localscratchOpt.strip().lower()
    if localscratchOpt=='true':       
        orig_list=copy.deepcopy(ext_list)
        try:
            with open(os.environ['PBS_NODEFILE'],'r') as nf:
                nf_string = nf.read()
            nf_list = list(set([x for x in nf_string.split('\n') if len(x.strip())!=0]))
	    try:
		    temp_dir= os.environ['TMPDIR']
	    except:
		    return

            work_dir=oneCalc['_AFLOWPI_FOLDER_']
#            ext_list=['wfc','hub','igk','restart','save','mix']

            prefix = AFLOWpi.retr.__prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])

            ext_list=['%s.%s' % (prefix,x) for x in ext_list]

            logging.info
        except Exception,e:
                  logging.debug(e)


        try:
            proc_by_node = [x for x in nf_string.split('\n') if len(x.strip())!=0]
            num_procs=len(proc_by_node)



            for ext in range(len(ext_list)):
                logging.info('copying %s files from:\nfrom : %s\nto   : %s' % (orig_list[ext],work_dir,temp_dir))
                proc_pool = multiprocessing.Pool(processes=num_procs)
                for proc in range(len(proc_by_node)):
                    try:

				    
			    try:
				    if orig_list[ext] in ['save','occup','update','bgfs']:
					    file_name=os.path.join(work_dir,'%s'%(ext_list[ext]))
				    elif orig_list[ext] in ['restart','restart_k','igk','restart_scf']:
					    if proc==0:
						    file_name=os.path.join(work_dir,'%s'%(ext_list[ext]))
					    else:
						    file_name=os.path.join(work_dir,'%s%s'%(ext_list[ext],proc+1))

				    else:
					    file_name=os.path.join(work_dir,'%s%s'%(ext_list[ext],proc+1))

			    except Exception,e:
				    continue

                            #add transfer to pool

			    """assuming the first file is in the folder we're executing this from"""
			    if not os.path.exists(file_name) and proc==0:
				    '''if we don't find the first skip to the nxt type'''
				    break

			    res = proc_pool.apply_async(AFLOWpi.prep.__scp_wrapper, (file_name,temp_dir,proc_by_node[proc])) 
			    if orig_list[ext] in ['save','occup','update','bgfs']:
				    break
		    except Exception,e:
			    logging.debug(e)
			    pass

                #closes pool for this specific ext type and joins to wait for all procs to finish
                proc_pool.close()
                proc_pool.join()


        except Exception,e:
            logging.debug(e)



def __get_tempdir():
    localscratchOpt=AFLOWpi.prep.__ConfigSectionMap('cluster','localscratch')
    localscratchOpt=localscratchOpt.strip().lower()
    if localscratchOpt=='true':        
        try:
            temp_dir= os.environ['TMPDIR']
            temp_dir= os.environ['TMPDIR']
        except Exception,e:
            logging.debug(e)
            temp_dir="./"
        if len(temp_dir.strip())==0:
            temp_dir="./"
    else:
        temp_dir="./"

    return temp_dir


def __setup_local_scratch(oneCalc,ID):
    localscratchOpt=AFLOWpi.prep.__ConfigSectionMap('cluster','localscratch')
    localscratchOpt=localscratchOpt.strip().lower()
    if localscratchOpt=='true':        
        temp_dir= AFLOWpi.prep.__get_tempdir()
        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','wfcdir',repr(temp_dir))
        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','outdir',repr(temp_dir))

#        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','wf_collect','.false.')
#        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','disk_io',"'medium'")

        
        #check if we're restarting and if so transfer the needed files to the local scratch for each respective node


#	AFLOWpi.prep.__to_local_scratch(oneCalc,ID,ext_list=['save'])
#        if AFLOWpi.prep.__check_restart(oneCalc,ID):
#	AFLOWpi.prep.__to_local_scratch(oneCalc,ID)


#        try:
#            AFLOWpi.prep.__to_local_scratch(oneCalc,ID)
#        except Exception,e:
#            AFLOWpi.run._fancy_error_log(e)
    else:        
        temp_dir='./'
        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','wfcdir',repr(temp_dir))
        oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'&control','outdir',repr(temp_dir))

    logging.info('setting wfc_dir to local scratch: %s' % temp_dir)

    return oneCalc,ID


#########################################################################################################################
#########################################################################################################################
## CALC SET FORMATION
#########################################################################################################################
#########################################################################################################################


import inspect
import itertools as it
def build_calcs(PARAM_VARS,build_type='product'):
        
	if build_type == 'zip' and len(PARAM_VARS) > 0:
             chkListLen = False
#             listLen = len(PARAM_VARS[0])
             listLen = max(len(item) for item in PARAM_VARS)
             
             for entry in PARAM_VARS:
                     if len(entry) != listLen and len(entry)>1:
                             chkListLen = False
                             break
                     else:
                             chkListLen = True


             if chkListLen ==True:
                 '''allow length of the entries in a zip to all be the''' 
                 '''same or len==1 and just repeat the len==1 items'''
                 for entry in range(len(PARAM_VARS)):
                     if len(PARAM_VARS[entry])==1:
                         PARAM_VARS[entry]=tuple(it.repeat(PARAM_VARS[entry][0], listLen))

                 logging.info('Creating calculations by zipping parameter lists')
                 print "Creating calculations by zipping parameter lists"
                 CALCS_ITER = it.izip(*PARAM_VARS)
             else:
#                     logging.error("Parameter lists are of different lengths, switching to product mode")
                     print "ERROR: Parameter lists are of different lengths."
                     logging.error("ERROR: Parameter lists are of different lengths on zip build mode. Exiting")
                     raise SystemExit
#                     CALCS_ITER = it.product(*PARAM_VARS)

        elif build_type == 'product' or len(PARAM_VARS) == 0:
             logging.info('Creating calculations in product mode')
             print "Creating calculations in product mode"
             CALCS_ITER = it.product(*PARAM_VARS)

        else:
            
             logging.info('ERROR: Unidentified build mode')
             print "ERROR: Unidentified build mode"
             raise SystemExit
	return CALCS_ITER
	


def scfs(aflowkeys,allAFLOWpiVars, refFile,pseudodir=None,build_type='product',build=True):
	"""
	Read a reference input file, and construct a set of calculations from the allAFLOWpiVars 
	dictionary defining values for the keywords in the reference input file. This will
	also create directory within the set directory for every calculation in the set.

	Arguments:
 	      allAFLOWpiVars (dict): a dictionary whose keys correspond to the keywords in the 
	                          reference input file and whose values will be used to 
   			          construct the set of calculations
  	      refFile (str): a filename as a string, a file object, or a string of the file
	                      that contains keywords to construct the inputs to the different
		              calculations in the set
	Keyword Arguments:
	      pseudodir (str): path of the directory that contains your Pseudopotential files
	                        The value in the AFLOWpi config file used will override this.
	      build_type (str): how to construct the calculation set from allAFLOWpiVars dictionary:
	                 
	                      zip | The first calculation takes the first entry from the list of 
			          | each of the keywords. The second calculation takes the second
			          | and so on. The keywords for all lists in allAFLOWpiVars must be
			          | the same length for this method.

			      product | Calculation set is formed via a "cartesian product" with 
			              | the values the list of each keyword combined. (i.e if 
			              | allAFLOWpiVars has one keyword with a list of 5 entires and
				      | another with 4 and a third with 10, there would be 2000
			 	      | calculations in the set formed from them via product mode.  

	      build (bool): <DEFUNCT OPTION. NEEDS REMOVAL>
	
	Returns:
              A dictionary of dictionaries containing the set of calculations.
			 
	"""
	filterFunction=None
        try:
            with open(refFile,'r') as refFileObj:
                refFile = refFileObj.read()
        except:
            try:
                refFile = refFile.read()
            except:
                refFile=refFile



        workdir = __ConfigSectionMap('prep','workdir')
        if os.path.isabs(workdir) == False:
            configFileLocation = __getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            workdir =  os.path.join(configFileLocation, workdir)


        PROJDIR=aflowkeys['project']
        PROJECT=PROJDIR
        SET=aflowkeys['set']
        if pseudodir==None:
            pseudodir= __ConfigSectionMap('prep','pseudodir')
            if os.path.isabs(pseudodir) == False:
                configFileLocation = __getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                pseudodir =  os.path.normpath(os.path.join(configFileLocation, pseudodir))
        else:
            if os.path.isabs(pseudodir) == False:
                pseudodir =  os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), pseudodir))

	inputfile = refFile
	logging.debug("Entering prepAllCalcs")
	PARAM_VARS = []
	PARAM_LABELS = []
	DICT = collections.OrderedDict()

	for k,v in allAFLOWpiVars.items():
		if v != None:
			PARAM_VARS.append(v)
			PARAM_LABELS.append(k)
			
		

	#converts the values in allvars input from floats or ints to strings for processing with re module
        for entry in range(len(PARAM_VARS)): 
		if len(re.findall(PARAM_LABELS[entry],inputfile))==0:
			print 'INVALID REF INPUT FILE: KEYWORD %s NOT FOUND. EXITING.' % PARAM_LABELS[entry]
			raise SystemExit
		if type(PARAM_VARS[entry][0]) != type(''):

			try:
				PARAM_VARS[entry]=tuple([str(subEntry) for subEntry in PARAM_VARS[entry]])
			#if the input can't be converted into a string something in the input for that variable is wrong and the framework force exits
			except Exception: 
				print "ERROR: Type of input for variable %s is invalid. Please check your input." % PARAM_LABELS[entry]
				raise SystemExit
        try:


	    CALCS_ITER = build_calcs(PARAM_VARS,build_type)
            CALCS_TUPLE = tuple(CALCS_ITER)

            ak = '_AFLOWPI_PSEUDO_DIR_'
            av = os.curdir
            DICT.update({ak:av})

            inputfile = re.sub(ak,av,inputfile)
            ak = '_AFLOWPI_OUTDIR_'
            av = os.curdir
            DICT.update({ak:av})
            inputfile = re.sub(ak,av,inputfile)

            hash_string = ''
            calcs = collections.OrderedDict()
            configFileName = __getConfigFile()
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


	if CALCS_TUPLE[0]:

                index=0
		for i in CALCS_TUPLE:
			inputfile2 = inputfile
			D = copy.deepcopy(DICT)
			for k,v in enumerate(i):
                            try:
				key = PARAM_LABELS[k]
				D[key] = v
				inputfile2 = re.sub(key,v,inputfile2)

                            except Exception,e:
                                AFLOWpi.run._fancy_error_log(e)
#########################################################################################################################		
                            if re.search(r'_AFLOWPI_[A-Z][0-9]*_', key):
                                
                                #make sure user isn't trying to use an equation for atomic species.
                                if len(re.findall(ur'===[^,\n]*===',v))!=0:
                                    print """%s is not an acceptable format for species name. Examples of acceptable format are: 'Ga',"O","Sn",'K'..etc""" % v
                                    logging.error("""%s is not an acceptable format for species name. Examples of acceptable format are: 'Ga',"O","Sn",'K'..etc""" % v)



                                ap = key[:-1]+'PSEUDO_'
                                try: 
                                        speciesRe = ''.join([i for i in v if not i.isdigit()])

                                        vp = __getPseudofilename(speciesRe,pseudodir)
                                        
                                        
                                except Exception,e:
                                        print 'Could not get Pseudopotential filename for %s from the directory %s' % (v,pseudodir)
                                        logging.error('Could not get Pseudopotential filename for %s from the directory %s' % (v,pseudodir))
                                
                                try:
                                    inputfile2 = re.sub(ap,vp,inputfile2)		

                                    am = key[:-1]+'MASS_'
                                    vm = __getAMass(speciesRe)

                                    inputfile2 = re.sub(am,str(vm),inputfile2)		
                                    

                                    D.update({am:vm,ap:vp})
                                    D.update({'_AFLOWPI_CONFIG_':configFileName})
                                except Exception,e:
                                    AFLOWpi.run._fancy_error_log(e)


#########################################################################################################################

                        try:

                            inputfile2 = __resolveEqualities(inputfile2)

                            inputfile2 = __cleanInputStringSCF(inputfile2)                            
                            calc_label = __hash64String(inputfile2)

                            kp = '_AFLOWPI_PREFIX_'
			    calc_label+='_01'
                            vp ='_'+calc_label
			    prefix=vp

#                            inputfile2 = re.sub(kp,vp,inputfile2)		
                            D.update({kp:vp})
                            index+=1
                            D.update({'_AFLOWPI_INDEX_':index})

                            D.update({'_AFLOWPI_FOLDER_':os.path.join(workdir,PROJECT,SET,'%s_%s_%04d' % (PROJDIR, SET, D['_AFLOWPI_INDEX_']))})
                            inputStringDict = AFLOWpi.retr.__splitInput(inputfile2)
                            inputStringDict['&control']['prefix']=repr(prefix)
                            inputfile2 = AFLOWpi.retr.__joinInput(inputStringDict)
                            D.update({'PROJECT':aflowkeys['project']})
                            D.update({'SET':aflowkeys['set']})
                            D.update({'__refFile__':refFile})
                            D.update({'_AFLOWPI_INPUT_':inputfile2})
			    D['__chain_index__']=1
#                            try:
#                                D['__chain_index__']+=1
#                            except:
#                                
                            if calc_label not in calcs.keys():
                                calcs[calc_label] = copy.deepcopy(D)

                            D.clear()

                        except Exception,e:
                            AFLOWpi.run._fancy_error_log(e)

	else:
            try:
		print 'Single calculation: do you need AFLOWpi?'
		logging.info('Single calculation: do you need AFLOWpi?')

                inputfile = __resolveEqualities(inputfile)

                try:
                    inputfile = __cleanInputStringSCF(inputfile)
                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)

		calc_label = __hash64String(inputfile)

		kp = '_AFLOWPI_PREFIX_'
		calc_label+='_01'
		vp ='_'+calc_label
#		inputfile = re.sub(kp,vp,inputfile)		
		prefix=vp
		inputStringDict = AFLOWpi.retr.__splitInput(inputfile)
		inputStringDict['&control']['prefix']=repr(prefix)
		inputfile = AFLOWpi.retr.__joinInput(inputStringDict)

		DICT.update({kp:vp})
                DICT.update({'_AFLOWPI_INDEX_':1})
                DICT.update({'_AFLOWPI_CONFIG_':__getConfigFile()})
                DICT.update({'_AFLOWPI_FOLDER_':os.path.join(workdir,PROJECT,SET,'%s_%s_%04d' % (PROJDIR, SET, DICT['_AFLOWPI_INDEX_']))})

		DICT['__chain_index__']=1
#                except:
#			DICT['__chain_index__']=0
		DICT.update({'_AFLOWPI_INPUT_':inputfile}) 

                if calc_label not in calcs.keys():
                    calcs[calc_label] = copy.deepcopy(DICT)

            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)



        for k,v in calcs.iteritems():
            calcs[k]['prev']= []
            calcs[k]['__calcVarName__']= 'oneCalc'
            calcs[k]['__IDVarName__']= 'ID'

        '''delete the "!" for the vacancies'''
        calcsCopy=copy.deepcopy(calcs)
        for k,v in calcsCopy.iteritems():
            oneCalcCopy=copy.deepcopy(v)
            for l,m in oneCalcCopy.iteritems():
                if re.search(r'_AFLOWPI_[A-Z]\d*_',l):
                    if m.strip()=='!':
                        del calcs[k][l]


        for ID,oneCalc in calcs.iteritems():
            in_dict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
	    in_dict['&control']['outdir']="'./'"
	    in_dict['&control']['wfcdir']="'./'"
	    oneCalc['_AFLOWPI_INPUT__']=AFLOWpi.retr.__joinInput(in_dict)
	    
            calcs[ID]['__execCounter__']=0
            calcs[ID]['__qsubFileName__']=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.qsub' % ID) 
            calcs[ID]['__execCounter__']=0
            calcs[ID]['__walltime_dict__']={'start':40000000000.0,'walltime':5000000}
            calcs[ID]['__execCounterBkgrd__']=0
            calcs[ID]['__status__']=collections.OrderedDict({"Start":False,'Complete':False,"Restart":0,"Error":'None'})

        

        '''put the initial set of calcs in a log'''
#        updatelogs(calcs,'step_01',runlocal=True)


        calcs = calcs_container(calcs)

	
        if build==True:
            maketree(calcs, pseudodir=pseudodir)
        else:
            pass
                

#	if run==True and build==True:
#		AFLOWpi.run.scf(calcs)

	return calcs


def extractvars(refFile):
	"""
	Read refFile and return an empty dictionary with all the keys and None values
	
	Arguments:
     	      refFile (str): filename of the reference input file for the frame

	Returns:
	      A dictionary containing the keyword extracted from the reference input file 
	      as keys with None for values..

	"""

        try:
            with open(refFile,'r') as refFileObj:
                refFile = refFileObj.read()
        except:
            try:
                refFile = refFile.read()
            except:
                refFile=refFile
        

	rex1 = re.compile('_AFLOWPI_.+_')
	lvar = rex1.findall(refFile)

	allAFLOWpiVars = {}
	for k in lvar:
		for kk in re.split(' +',k):
			allAFLOWpiVars.update({kk:None})
	logging.debug('Exiting extractAFLOWpiVars')
# we should counsider adding consistency tests to checl that the ref file is good and working.

	return allAFLOWpiVars




def calcFromFile(aflowkeys,fileList,reffile=None,pseudodir=None,build=True,workdir=None,keep_name=False,clean_input=True):
	"""
	Reads in a string of an QE input file path, a string of an QE input, a file object of a 
	QE input or a list of them and attempts to fill create a calculation from them. If they
	are missing things such as k_points card, they are automtically generated. 


	Arguments:
	      aflowkeys (dict): a dictionary generated by AFLOWpi.prep.init
	      fileList (list): a string of an QE input file path, a string of an QE input, a file 
    	                        object of a QE input or a list of them
        Keyword Arguments:
	      reffile (str): a partially filled QE input file used in case the input(s) in fileList
	                      are missing. i.e. wfc cutoff. If the names of the Pseudopotential files 
			      are not included in the input(s) in fileList, they are chosen depending
			      on the pseudodir chosen and included when the calculation set is formed.
	      workdir (str): a string of the workdir path that be used to override what is in the 
	                      config file used when initating the AFLOWpi session
	      pseudodir (str): a string of the pseudodir path that be used to override what is in 
	                        the config file used when initating the AFLOWpi session
	      build (bool): <DEFUNCT OPTION. NEEDS REMOVAL>

	"""

        returnDict=collections.OrderedDict()
        index=0

        if type(fileList)==type('aString'):
            fileList=[fileList]

        for inputFile in fileList:
            try:
                holder=inputFile
                with open(inputFile,'r') as inputFileObj:
                    inputFile = inputFileObj.read()
		file_name=os.path.basename(holder)
		temp_file_hash_name='.'.join(file_name.split('.')[:-1])
            except:
                try:
                    inputFile = inputFile.read()
                except:
                    inputFile=inputFile

            '''try to open the reffile '''
            if reffile!=None:
                if type(reffile)==type('strang'):
			refFileStr=reffile
		else:
			try:
                    
				with open(reffile,'r') as inputFileObj:
					refFileStr = inputFileObj.read()
					refFileStr = __removeComments(refFileStr)
					'''get a dict of the tokenized ref file to fill in the parts missing from the files in filelist'''
			except Exception,e:
				print e
				'''if it doesn't work just return a blank dict for the ref file'''
				refDict=collections.OrderedDict() 
		try:
		       refDict = AFLOWpi.retr.__splitInput(refFileStr)
		except:
		       refDict=collections.OrderedDict() 

	    else:
		    '''if the user doesn't input a ref file then just return a blank dict'''
		    refDict=collections.OrderedDict() 
			

                    
                    
            inputFile = __removeComments(inputFile)

            if workdir==None:            
                workdir = __ConfigSectionMap('prep','workdir')

            if os.path.isabs(workdir) == False:
                configFileLocation = __getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                workdir =  os.path.join(configFileLocation, workdir)



    #        configFileName = os.path.dirname(configFileLocation)
            PROJDIR=aflowkeys['project']
            PROJECT=PROJDIR
            SET=aflowkeys['set']
            if pseudodir==None:
                pseudodir= __ConfigSectionMap('prep','pseudodir')
                if os.path.isabs(pseudodir) == False:
                    configFileLocation = __getConfigFile()
                    configFileLocation = os.path.dirname(configFileLocation)
                    pseudodir =  os.path.normpath(os.path.join(configFileLocation, pseudodir))
            else:
                if os.path.isabs(pseudodir) == False:
                    pseudodir =  os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), pseudodir))

            inputfile = inputFile
            configFileName = __getConfigFile()

            DICT = collections.OrderedDict()
            counter=1
            elementList = list(string.ascii_uppercase)

            numOfEach = AFLOWpi.retr.__getAtomNum(inputfile,strip=False)
            numOfEachStripped = AFLOWpi.retr.__getAtomNum(inputfile,strip=True)
#	    print numOfEachStripped
	    num_map={}
	    map_counter=1
	    for species,num in numOfEachStripped.iteritems():
		    num_map[species]=map_counter
		    map_counter+=1
	    
            atomicSpeciesList=[]
	    element_list=list(string.ascii_uppercase)
            for species,num in numOfEach.iteritems():
                try:
                    strippedSpecies = species.strip("0123456789")
                    numSpecies = species.strip(string.ascii_uppercase+string.ascii_lowercase)
#
		    counterIndex=counter/26
#                    if counter/26==0:
#                        counterString=''

			
#                    else:
#                        counterString=str((counter)%26) 
		    AFLOWpi_name = element_list[num_map[strippedSpecies]]+numSpecies
                    aSpec='_AFLOWPI_%s_' % AFLOWpi_name

                    vSpec=strippedSpecies
                    vp = __getPseudofilename(strippedSpecies,pseudodir)            
                    ap = '_AFLOWPI_%sPSEUDO_' % AFLOWpi_name
                    am = '_AFLOWPI_%sMASS_' % AFLOWpi_name
                    vm = __getAMass(strippedSpecies)        
                    atomicSpeciesList.append('%s %s %s' % (species,vm,vp))
                    
                    DICT.update({am:vm,ap:vp,aSpec:species})
                    counter+=1

                except Exception,e:
			print e
			raise SystemExit
#                    print 'Could not get Pseudopotential filename for %s from the directory %s' % (species,pseudodir)
 #                   logging.error('Could not get Pseudopotential filename for %s from the directory %s' % (species,pseudodir))

            atomicSpecString='\n'.join(atomicSpeciesList)+'\n'

            inputCalc = AFLOWpi.retr.__splitInput(inputfile)

#            for namelist in refDict.keys():
#                for param,value in namelist.iteritems():
                    


            for namelist,entries in refDict.iteritems():
                try:
                    if namelist not in inputCalc.keys():
                        '''if an input namelist from the ref isn't in the input file add it '''
                        inputCalc[namelist]=collections.OrderedDict()
                    for variable,value in entries.iteritems():
                        if variable not in inputCalc[namelist].keys():
                            '''if a variable from the namelist in the ref isn't in the input file add it'''
                            inputCalc[namelist][variable]=value
                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)
                    pass

            if 'ATOMIC_SPECIES' not in inputCalc.keys() or '__content__' not in inputCalc['ATOMIC_SPECIES'].keys() or inputCalc['ATOMIC_SPECIES']['__content__'].strip()=='':
                inputCalc['ATOMIC_SPECIES']=collections.OrderedDict()
                inputCalc['ATOMIC_SPECIES']['__modifier__']=''

                inputCalc['ATOMIC_SPECIES']['__content__']=atomicSpecString



            if '&control' not in inputCalc.keys():
                inputCalc['&control']=collections.OrderedDict()
#                inputCalc['&control']['calculation']='"scf"'
#                inputCalc['&control']['tstress']='.TRUE.'
#                inputCalc['&control']['tprnfor']='.TRUE.'
            inputCalc['&control']['PSEUDO_DIR']="'./'"
            inputCalc['&control']['outdir']="'./'"
	    inputCalc['&control']['wfcdir']="'./'"

            inputCalc['&control']['restart_mode']="'from_scratch'"
	    
            if '&electrons' not in inputCalc.keys():
                inputCalc['&electrons']=collections.OrderedDict()
                inputCalc['&electrons']['diagonalization']='"david"'
                inputCalc['&electrons']['mixing_mode']='"plain"'
                inputCalc['&electrons']['mixing_beta']='0.7'
                inputCalc['&electrons']['conv_thr']='1.0d-8'

                if 'ecutwfc' not in inputCalc['&system'].keys():
                    try:
                        
                        ecutwfc = inputCalc['&system']['ecutwfc']
                        ecutrho=str(int(float(ecutwfc)*4))
                        inputCalc['&system']['ecutrho']=ecutrho
                    except: 
                        ecutwfc='200'
                        ecutrho='800'
                
                        DICT.update({'_AFLOWPI_ECUTW_':ecutwfc})
                        DICT.update({'_AFLOWPI_ECUTR_':ecutrho})





            inputCalc = AFLOWpi.retr.__orderSplitInput(inputCalc)
            inputfile = AFLOWpi.retr.__joinInput(inputCalc)


            DICT.update({'_AFLOWPI_CONFIG_':configFileName})

            ak = '_AFLOWPI_PSEUDO_DIR_'
            av = './'
            DICT.update({ak:av})
            ak = '_AFLOWPI_OUTDIR_'
            av = './'
            DICT.update({ak:av})
            index+=1


            try:
                if '__content__' not in inputCalc['K_POINTS'].keys() or inputCalc['K_POINTS']['__content__'].strip()=='':
                        '''if we don't have a kpoints card make a generic one based on lattice vecs and MP Grid'''

                        cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(inputfile)

                        inputCalc = AFLOWpi.retr.__splitInput(inputfile)


                        kpointString=getMPGrid(cellParamMatrix,offset=True)
			if kpointString=='':
				inputCalc['K_POINTS']['__modifier__']='{gamma}'
				inputCalc['K_POINTS']['__content__']=''
			else:
				inputCalc['K_POINTS']=collections.OrderedDict()
				inputCalc['K_POINTS']['__modifier__']='{automatic}'
				inputCalc['K_POINTS']['__content__']=kpointString

                        inputfile = AFLOWpi.retr.__joinInput(inputCalc)

            except Exception,e:
                print e
                AFLOWpi.run._fancy_error_log(e)
                

	    if clean_input==True:
		    inputfile = __cleanInputStringSCF(inputfile)                            



            '''make the ID'''
            try:
#                '''try to see if there's a global counter already in play'''
#		try:
#			chain_ind = int(globals()['__TEMP__INDEX__COUNTER__'])
#		except:
		chain_ind=1


		ext ='_%02d'%chain_ind

            except Exception,e:
                print e
                ext=''
            


            calc_label = __hash64String(inputfile)+ext
	    if keep_name==True:
		    try:
			    calc_label=temp_file_hash_name+ext
		    except:
			    pass
            kp = '_AFLOWPI_PREFIX_'

            vp ='_'+calc_label
            DICT.update({kp:vp})

            inputCalc = AFLOWpi.retr.__splitInput(inputfile)

            inputCalc['&control']['prefix']=repr(vp)


            inputfile = AFLOWpi.retr.__joinInput(inputCalc)


            DICT['prev']=[]

            DICT['__execCounter__']=0
            DICT['__execCounterBkgrd__']=0


            
            DICT.update({'_AFLOWPI_INDEX_':index})
	    DICT.update({'__chain_index__':1})

	    if clean_input==True:
		    inputfile = __cleanInputStringSCF(inputfile)                            


            DICT.update({'_AFLOWPI_INPUT_':inputfile})
            DICT.update({'__refFile__':inputfile})
            DICT.update({'_AFLOWPI_FOLDER_':os.path.join(workdir,PROJECT,SET,'%s_%s_%04d' % (PROJDIR, SET, DICT['_AFLOWPI_INDEX_']))})
            DICT['__qsubFileName__']=os.path.join(DICT['_AFLOWPI_FOLDER_'],'_%s.qsub' % calc_label) 
            DICT['__calcVarName__']= 'oneCalc'
            DICT['PROJECT']= PROJECT
            DICT['SET']= SET
            DICT['__IDVarName__']= 'ID'
            DICT['__status__']=collections.OrderedDict({"Start":False,'Complete':False,"Restart":0,"Error":'None'})
            dictCopy=copy.deepcopy(DICT)
            if calc_label not in returnDict:
                returnDict[calc_label]=dictCopy

        


        if build==True:
            maketree(returnDict, pseudodir=pseudodir,workdir=workdir)

        else:
            pass

        for ID,oneCalc in returnDict.iteritems():
            __saveOneCalc(oneCalc,ID)
        return returnDict



def getMPGrid(primLatVec,offset=True,string=True):
    try:
        '''just starting on this. going to refine later'''

        kpointList=[]
        a,b,c,alpha,beta,gamma =  AFLOWpi.retr.free2abc(primLatVec,cosine=False,bohr=False,string=False)
        kpointList.append(int(numpy.floor(200.0/(a*2*numpy.pi))))
        kpointList.append(int(numpy.floor(200.0/(b*2*numpy.pi))))
        kpointList.append(int(numpy.floor(200.0/(c*2*numpy.pi))))



        if offset==True:
            kpointList.extend([1,1,1,])
        else:
            kpointList.extend([0,0,0,])


	if kpointList[:3]==[1,1,1]:
		return ''
	

        if string==True:
            return ' '.join([str(x) for x in kpointList])
        else:
            return kpointList
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        return '2 2 2 1 1 1'
        
#########################################################################################################################
#########################################################################################################################
## CONFIG UTILS
#########################################################################################################################
######################################################################################################################### 


def __getConfigFile():
    return __config__

def __forceGlobalConfigFile(configFile):
    global __config__
    __config__=configFile



def __ConfigSectionMap(section,option,configFile=None):
    index=1
#    try:
#	    index= globals()['__TEMP__INDEX__COUNTER__']
#    except:
#	    pass

    Config = ConfigParser.ConfigParser()
    config = __getConfigFile()
    if configFile == None:
        Config.read(config)
    else:
        Config.read(configFile)        

    dict1 = {}

    try:
        options = Config.options(section)
    except ConfigParser.NoSectionError:
        logging.info('No section in %s called: %s' % (config,section))
        return ''



    if option in options:
        try:
            returned_option = Config.get(section, option)


        except:
            logging.info('No section in %s in section %s called: %s' % (config,section,option))
            returned_option=''
    else:
            returned_option=''
    try:
        section='step_%02d'%int(index)
        returned_option = Config.get(section, option)
    except Exception,e:
        pass

    return returned_option


#########################################################################################################################
#########################################################################################################################
## INTERFACE
#########################################################################################################################
######################################################################################################################### 


def ConfigSectionMap(section,option,configFile=None):
	return __ConfigSectionMap(section,option,configFile=configFile)
#########################################################################################################################
class init:

        def __init__(self,PROJECT, SET='', AUTHOR='', CORRESPONDING='', SPONSOR='',config='',workdir=None):
		self.project=PROJECT
		self.set=SET
		self.config=config
		self.author=AUTHOR
		self.corresponding=CORRESPONDING
		self.sponsor=SPONSOR

		self.workdir=workdir	
#		if '__global__config__flag__' in globals():
#			print 'Only one session per script. Exiting'
#			raise SystemExit
		globals()['__global__config__flag__']=True

		self.keys = AFLOWpi.prep.init__(PROJECT, SET=SET, AUTHOR=AUTHOR, CORRESPONDING=CORRESPONDING, SPONSOR=SPONSOR,config=config,workdir=workdir)



        def items(self):
                return self.keys.items()

        def iteritems(self):
                return ((i,j) for i,j in self.keys.iteritems())

        def keys(self):
                return self.keys.keys()

        def values(self):
                return self.keys.values()

        def __getitem__(self,index):
                return self.keys[index]

        def __delitem__(self,index):
                del self.keys[index]

        def __str__(self):
                return repr(self.keys)

        def __repr__(self):
                return self.__str__()

        def __len__(self):
                return len(self.keys)

	def scfs(self,allAFLOWpiVars, refFile,name='first',pseudodir=None,build_type='product',build=True,run=True):
		"""
		A wrapper method to call AFLOWpi.prep.scfs to form the calculation set. This will
		also create directory within the set directory for every calculation in the set.

		Arguments:
		      allAFLOWpiVars (dict): a dictionary whose keys correspond to the keywords in the 
				          reference input file and whose values will be used to 
				          construct the set of calculations
		      refFile (str): a filename as a string, a file object, or a string of the file
				      that contains keywords to construct the inputs to the different
				      calculations in the set
		Keyword Arguments:
		      pseudodir (str): path of the directory that contains your Pseudopotential files
				        The value in the AFLOWpi config file used will override this.
		      build_type (str): how to construct the calculation set from allAFLOWpiVars dictionary:

				      zip | The first calculation takes the first entry from the list of 
					  | each of the keywords. The second calculation takes the second
					  | and so on. The keywords for all lists in allAFLOWpiVars must be
					  | the same length for this method.

				      product | Calculation set is formed via a "cartesian product" with 
					      | the values the list of each keyword combined. (i.e if 
					      | allAFLOWpiVars has one keyword with a list of 5 entires and
					      | another with 4 and a third with 10, there would be 2000
					      | calculations in the set formed from them via product mode.  

		      build (bool): <DEFUNCT OPTION. NEEDS REMOVAL>
		      run (bool): <DEFUNCT OPTION. NEEDS REMOVAL>

		Returns:
		      A dictionary of dictionaries containing the set of calculations.

		"""

		scfs = AFLOWpi.prep.scfs(self.keys,allAFLOWpiVars, refFile,pseudodir=pseudodir,build_type=build_type,build=build)
		return scfs


	def from_file(self,fileList,reffile=None,pseudodir=None,build=True,workdir=None):
		"""
		Reads in a string of an QE input file path, a string of an QE input, a file object of a 
		QE input or a list of them and attempts to fill create a calculation from them. If they
		are missing things such as k_points card, they are automtically generated. 

		Arguments:
		      aflowkeys (dict): a dictionary generated by AFLOWpi.prep.init
		      fileList (str): a string of an QE input file path, a string of an QE input, a file 
				       object of a QE input or a list of them
		Keyword Arguments:
		      reffile (str): a partially filled QE input file used in case the input(s) in fileList
		                      are missing. i.e. wfc cutoff. If the names of the Pseudopotential files 
              			      are not included in the input(s) in fileList, they are chosen depending
				      on the pseudodir chosen and included when the calculation set is formed.
		      workdir (str): a string of the workdir path that be used to override what is in the 
				   config file used when initating the AFLOWpi session
		      pseudodir (str): a string of the pseudodir path that be used to override what is in 
				     the config file used when initating the AFLOWpi session
		      build (bool): <DEFUNCT OPTION. NEEDS REMOVAL>

		"""

		scfs=AFLOWpi.prep.calcFromFile(self.keys,fileList,reffile=reffile,pseudodir=pseudodir,build=build,workdir=workdir)
		return calcs_container(scfs)

	def load(self,step=1):
		"""
		Loads the calc logs from a given step

		Arguments:
		      step (int): the step of the calculation for whose calclogs are to be loaded
		
		Returns:
		      calcs (dict): the loaded calc logs

		"""
		loaded_log_name='step_%02d'%step
		loaded_calcs = AFLOWpi.prep.loadlogs(PROJECT=self.project,SET=self.set,logname=loaded_log_name,config=self.config)
		return calcs_container(loaded_calcs)

	def status(self,status={},step=0,negate_status=False):
		"""
		Loads the calc logs from a given step

		Arguments:
		      step (int): The step of the calculation for whose calclogs are to be loaded.
              		          If no step is specified then it will default to load calculations
				  from all steps with the chosen status.

		Keyword Arguments:
		      status (dict): key,value pairs for status type and their value to filter on.
		                     i.e. status={'Finished':False} 
		      negate_status (bool): filter on the opposite of the status filters
		
		Returns:
		      calcs (dict): the loaded calcs for one or more steps with the given status

		"""

#		loaded_log_name='step_%02d'%step
		loaded_calcs = AFLOWpi.retr.checkStatus(self.project,SET=self.set,config=self.config,step=step,status=status,negate_status=negate_status)
#		loaded_calcs = AFLOWpi.prep.loadlogs(PROJECT=self.project,SET=self.set,logname=loaded_log_name,config=self.config)
		full_set=collections.OrderedDict()
		for i in loaded_calcs:
			full_set.update(i)
		
		return calcs_container(full_set)


		



def __getLoglevel():
    """
    Checks config file for loggings level. 
    If none specified it defaults to logging.INFO

    Arguments:
          None

    Returns:
          loglevel (logging.loglevel): loglevel object associated with the logging
	                               level string found in the config file

    """
    

    loglevel = __ConfigSectionMap('prep','loglevel').upper()
    if loglevel=='DEBUG':
        loglevel=logging.DEBUG
    elif loglevel == 'WARNING':
        loglevel=logging.WARNING
    elif loglevel == 'ERROR':
        loglevel=logging.ERROR
    elif loglevel == 'CRITICAL':
        loglevel=logging.CRITICAL
    else:
        loglevel=logging.INFO

    return loglevel

#########################################################################################################################



def init__(PROJECT, SET='', AUTHOR='', CORRESPONDING='', SPONSOR='',config='',workdir=None):

	"""
	Initializes the frame
	
	Arguments:
	 PROJECT (str): Name of project
         SET (str): Name of set 
	 author (str): Name of author
	 CORRESPONDING (str): Name of corresponding 
	 SPONSOR (str): Name of sponsor
	 
	e.g. initFrame('LNTYPE','', 'MF', 'marco.fornari@cmich.edu','DOD-MURI'). Return the AFLOKEYS dictionary.
	"""  
        configFile = os.path.abspath(os.path.join(os.curdir,config))
        print 'Reading from Config: %s' % configFile

        try:
            globals()['__INITIAL_EXECPATH__']=os.path.realpath(__file__)
        except NameError,e:
            globals()['__INITIAL__EXECPATH__']=os.path.dirname(sys.argv[0])

	AFLOWpi.prep.__forceGlobalConfigFile(configFile)

# for classification purpose (AFLOWLIB like): 
#	mkdir PROJECT, the bottom of the directory tree
#	define and export a dictionary: author, corresponding, sponsor, date
#	define a (mkdir) SET Layer, second layer of the directory tree, may or may not exist physically
# start the logging
        if workdir==None:
            workdir = __ConfigSectionMap('prep','workdir')
        if os.path.isabs(workdir) == False:
                    configFileLocation = __getConfigFile()
                    configFileLocation = os.path.dirname(configFileLocation)
                    workdir =  os.path.join(configFileLocation, workdir)


        if not os.path.exists(workdir):
            try:
                os.makedirs(workdir)
            except Exception,e:
                logging.warning('work directory %s does not exist and attempting to create it failed. Exiting AFLOWpi' % workdir)
                SystemExit
        if workdir==None:
            workdir=os.curdir
        
	path = os.path.join(workdir,PROJECT,SET)

#	os.chdir(os.path.abspath(workdir))
	AFLOWpidir = os.path.join(path,'AFLOWpi')

	if os.path.exists(AFLOWpidir):
                logging.info("Project exists: %s" %  path)
		print "Project exists: ", path
 		pass
	else:
            try:
                try:
			os.makedirs(path)
		except:
			pass
		AFLOWpidir = os.path.join(path,'AFLOWpi')
		os.makedirs(AFLOWpidir)
		calclogDir=os.path.join(AFLOWpidir,'calclogs')
		os.makedirs(calclogDir)
            except:
                pass
	if os.path.exists(AFLOWpidir):
		logfile = os.path.join(AFLOWpidir,'LOG.log')
			

		logging.basicConfig(filename=logfile,format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',level=__getLoglevel())
#		logging.setLoggerClass(AFLOWpi.prep.colorful_logging)
		logging.info('Start logging for AFLOWpi for project %s in %s', PROJECT, AFLOWpidir )

                loglevel = ''
                if __ConfigSectionMap('prep','loglevel') == '':
                    loglevel = 'INFO'
                else:
                    loglevel = __ConfigSectionMap('prep','loglevel').upper()
		logging.info('Log Level set at: logging.%s' % loglevel)
	
	date = datetime.date.isoformat(datetime.date.today())


	__aflowkeys__ = {}
	__aflowkeys__.update({'project':PROJECT, 'set':SET, 'author':AUTHOR, 'corresponding':CORRESPONDING, 'sponsor':SPONSOR})
	__aflowkeys__.update({'date':date})

        dest = os.path.join(AFLOWpidir,'CONFIG.config')        
        AFLOWpi.prep.__copyConfig(config,dest)
	


	AFLOWpi.prep.__forceGlobalConfigFile(dest)
	
	return __aflowkeys__

def __copyConfig(config,dest):
    '''
    Copies the config file used to the AFLOWpi directory when AFLOWpi 
    is initiated and resolves relative paths in the config file.

    Arguments:
          config (str): location of the config file to be copied to "AFLOWpi" directory


    '''



    try:
        if not os.path.isabs(config):

            execDirName =  os.path.dirname(os.path.realpath(__main__.__file__))
            config  =  os.path.normpath(os.path.join(execDirName,config))


    except Exception,e:
        print e

    configParse  = ConfigParser.ConfigParser()
    configParse.read(config)
    sectionList = configParse.sections()

    configFileDir = os.path.dirname(config)

    
    for section in sectionList:
        optionList = configParse.options(section)        
        for option in optionList:
            if option == '' or option == 'paodir' or option == 'pseudodir' or option == 'espressodir' or option == 'wantdir' or option == 'workdir' or option == 'enginedir' or option == 'jobreffile':
                possiblePath = configParse.get(section, option)
                possiblePath  =  os.path.normpath(os.path.join(configFileDir, possiblePath))       

                if os.path.abspath(possiblePath):
                    pass
                else:
                    try:
                        possiblePath =  os.path.normpath(os.path.join(configFileDir,possiblePath))            
                    except Exception,e:
                        print e
                if os.path.isfile(possiblePath) or os.path.isdir(possiblePath):
                    pass
                else:
                    print 'entry in %s in section %s does not exist. check your config file. If you do not need/have this directory then either comment out or delete this entry in your config and run your script again.' % (option,section)
                    logging.error('entry in %s in section %s does not exist. check your config file. If you do not need/have this directory then either comment out or delete this entry in your config and run your script again.' % (option,section))
                    raise SystemExit
                configParse.set(section, option, possiblePath)



                

    with open(dest,'w') as configDest:
        configParse.write(configDest)



class calcs_container:
        def __init__(self,dictionary):

                self.int_dict=dictionary
		self.type='scf'
		self.step_index=0
		self.load_index=0
		self.plot=plotter(self.int_dict)
		self.split_index=0
		self.initial_inputs=self.__getInitInputs()
		self.initial_calcs=[]
		self.split=True
		self.tight_banding=False
		self.scf_complete=False
		self.workflow=[]
		self.submit_flag=False
#	def implicit(self):



	def increase_step(func):
		@functools.wraps(func)
		def wrapper(*args,**kwargs):
			self.step_index+=1
			print self.step_index
			return func(*args, **kwargs)


        def items(self):
                return self.int_dict.items()

        def iteritems(self):
                return ((i,j) for i,j in self.int_dict.iteritems())

        def keys(self):
                return self.int_dict.keys()

        def values(self):
                return self.int_dict.values()

        def __getitem__(self,index):
                return self.int_dict[index]

        def __delitem__(self,index):
                del self.int_dict[index]

        def __str__(self):
                return repr(self.int_dict)

        def __repr__(self):
                return self.__str__()

        def __len__(self):
                return len(self.int_dict)

	def __getInitInputs(self):
		inputDict=collections.OrderedDict()
		for ID,oneCalc in self.int_dict.iteritems():
			inputDict[ID.split('_')[0]]=oneCalc['_AFLOWPI_INPUT_']

		return inputDict

	
	def conventional_cell_input(self):
		for ID,oneCalc in self.int_dict.iteritems():
			conv_input = AFLOWpi.retr.transform_input_conv(oneCalc,ID)
			self.int_dict[ID]['_AFLOWPI_INPUT_']=conv_input
			oneCalc['_AFLOWPI_INPUT_']=conv_input
			self.__saveOneCalc(oneCalc,ID)		     

			with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as newIn:
				newIn.write(conv_input)		


#	def submit(self):
#		AFLOWpi.run.submit()
#		AFLOWpi.run.self.initial_inputs[0]
	
	

	def get_initial_inputs(self):
		temp_dict=collections.OrderedDict()
		temp_dict=copy.deepcopy(self.int_dict)
		for ID,oneCalc in self.int_dict.iteritems():
			temp_dict[ID]['_AFLOWPI_INPUT_']=self.initial_inputs[ID.split('_')[0]]
		self.int_dict=copy.deepcopy(temp_dict)
		return self

	def __saveOneCalc(self,oneCalc,ID):
		with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.oneCalc' % ID),'wb') as oneCalcPickleFile:
			oldOneCalc = cPickle.dump(oneCalc,oneCalcPickleFile)    


	def resubmit(self,reset=True):
		AFLOWpi.run.reset_logs(self.int_dict)
		AFLOWpi.run.resubmit(self.int_dict)


	def split_chain(self,update_positions=True,update_structure=True):
		print '\nSPLITTING CHAIN\n'
		new_calcs=calcs_container(self.int_dict)	
		split_index=1		
		try:
			
			split_index=globals()['__SPLIT_INDEX__']+1
			globals()['__SPLIT_INDEX__']=split_index
		except Exception,e:
			globals()['__SPLIT_INDEX__']=split_index

		try:
			last_calcs = self.initial_calcs[-1]
		except:
			print 'you cannot split a chain with zero steps..exiting'
			raise SystemExit
		
		new_calcs.int_dict=last_calcs
		new_calcs.split=True
#		new_calcs.new_step(update_positions=True,update_structure=True,new_job=False)




		


#=======
#		new_calcs.new_step(update_positions=True,update_structure=True,new_job=False)	
#		pre = 'split_%02d' % split_index	
#		new_calcs.int_dict = AFLOWpi.prep.modifyInputPrefixPW(new_calcs.int_dict,repr(pre))
#>>>>>>> WanT_conduction

		return new_calcs


#	@increase_step
        def scf(self):
		self.scf_complete=True
		self.tight_banding==False
		self.load_index+=1
		self.type='scf'
		self.new_step(update_positions=True,update_structure=True)
		self.initial_calcs.append(self.int_dict)
		self.change_input('&control','calculation','"scf"')#,change_initial=False)
		calc_type='scf'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)

		AFLOWpi.run.scf(self.int_dict)	

#	@increase_step
	def relax(self):
		self.scf_complete=True
		self.tight_banding==False
#		self.step_index+=1
		self.type='relax'
		self.load_index+=1
		self.new_step(update_positions=True,update_structure=True,)
		self.initial_calcs.append(self.int_dict)
		self.change_input('&control','calculation','"relax"')#,change_initial=False)

		calc_type='Ionic Relaxation'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)

		AFLOWpi.run.scf(self.int_dict)	

#	@increase_step		
	def vcrelax(self):
		self.scf_complete=True
		self.tight_banding==False
#		self.step_index+=1
		self.type='vcrelax'
		self.load_index+=1
		self.new_step(update_positions=True,update_structure=True)
		self.initial_calcs.append(self.int_dict)
		self.change_input('&control','calculation','"vc-relax"')#,change_initial=False)

		calc_type='Variable Cell Relaxation'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)

		AFLOWpi.run.scf(self.int_dict)	
		
	def addToAll(self,block=None,addition=None):
		if block==None or addition==None:
			return
		for ID,oneCalc in self.int_dict.iteritems():
			AFLOWpi.prep.addToBlockWrapper(oneCalc,ID,block,addition)    
				
	
        def new_step(self,update_positions=True,update_structure=True,new_job=True,ext=''):
		self.step_index+=1
#		globals()['__TEMP__INDEX__COUNTER__']=self.step_index-1
		try:
			last_calcs = self.initial_calcs[-1]
		except:
			last_calcs=self.int_dict

		'add to list of steps on oneCalcs'


		self.workflow.append(self.type)

			
		self.int_dict = writeToScript('',last_calcs,from_step=self.step_index)
		'chooses if the positions and/or structs to update'
		if self.split==True:
			pre = 'split_%02d' % self.step_index
			self.int_dict = AFLOWpi.prep.modifyInputPrefixPW(self.int_dict,repr(pre))

			self.split=False

		self.int_dict=updateStructs(self.int_dict,update_positions=update_positions,update_structure=update_structure)

				


#		new_calcs = self.get_initial_inputs()
#		new_calcs.init_inputs=self.__()		

#		for k,v in self.int_dict.iteritems():
#			print v['_AFLOWPI_INDEX_']
#		except:
#			new_calcs.step_index=1

#		globals()['__TEMP__INDEX__COUNTER__']=new_calcs.step_index
		self.plot=plotter(self.int_dict)
#		return  new_calcs

	def change_pseudos(self,directory):
		if not os.path.isabs(directory):
			directory=os.path.join(os.path.dirname(os.path.realpath(__main__.__file__)),directory)

		
		self.int_dict = changeCalcs(self.int_dict,keyword='pseudos',value=directory)

#		self.intial_calcs = AFLOWpi.prep.changeCalcs(oneCalc,ID,keyword='calculation',value='scf')

	def change_input(self,namelist=None,parameter=None,value=None):
		self.int_dict=modifyNamelistPW(self.int_dict,namelist,parameter,value)

#		self.intitial_calcs = AFLOWpi.prep.modifyNamelistPW(self.initial_calcs,namelist,parameter,value)

#	@increase_step
	def phonon(self,nrx1=2,nrx2=2,nrx3=2,innx=2,de=0.01,mult_jobs=False,raman=False,LOTO=False,disp_sym=True,atom_sym=True,field_strength=0.01):

		self.tight_banding=False
		self.type='phonon'
		self.new_step(update_positions=True,update_structure=True)
		
		loadModString = '''AFLOWpi.run.prep_fd(__submitNodeName__,oneCalc,ID,nrx1=%i,nrx2=%i,nrx3=%i,innx=%i,de=%f,atom_sym=%s,disp_sym=%s)'''%(nrx1,nrx2,nrx3,innx,de,atom_sym,disp_sym)
		self.addToAll(block='PREPROCESSING',addition=loadModString)		

		force_pull_string="""AFLOWpi.prep.__addToAll(calc_subset,'RUN','AFLOWpi.run.__pull_forces(oneCalc,ID)')"""
		if raman==True:
			LOTO=True
			task_list=['AFLOWpi.run.scf(calc_subset)','AFLOWpi.run.__setup_raman(calc_subset,oneCalc,field_strength=%s)'%field_strength,force_pull_string,]
		elif LOTO==True:
			task_list=['AFLOWpi.run.scf(calc_subset)','AFLOWpi.run.__setup_raman(calc_subset,oneCalc,for_type="born",field_strength=%s)'%field_strength,force_pull_string,]
		else:
			task_list=['AFLOWpi.run.scf(calc_subset)',force_pull_string,]

		self.int_dict=AFLOWpi.prep.prep_split_step(self.int_dict,'AFLOWpi.prep.__one_phon(oneCalc,ID,subset_name="FD_PHONON")',subset_tasks=task_list,mult_jobs=mult_jobs,substep_name='FD_PHONON',keep_file_names=True,clean_input=True)
#		self.addToAll(block='RUN',addition='AFLOWpi.run.raman(oneCalc,ID)')				



		loadModString = '''AFLOWpi.run.__pp_phonon(__submitNodeName__,oneCalc,ID,de=%s,raman=%s,LOTO=%s,field_strength=%s)'''%(de,raman,LOTO,field_strength)
		self.addToAll(block='POSTPROCESSING',addition=loadModString)		

		calc_type='Phonon'
		if raman==True:
			calc_type+=' with Raman scattering'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)



#pyt	@increase_step(len(temperature))
	def transport(self,temperature=[300,],epsilon=True,run_bands=True):
		'''
		Wrapper method to call AFLOWpi.scfuj.prep_transport and AFLOWpi.scfuj.run_transport 
		in the high level user interface. Adds a new step to the workflow.
		
		

		Arguments:
		      self: the _calcs_container object

		Keyword Arguments:
		      epsilon (bool): if True episilon tensor will be computed 
		      temperature (list): list of temperature(s) at which to calculate transport properties

		Returns:
		      None

		'''		

		self.type='transport'



		self.new_step(update_positions=True,update_structure=True)
		self.initial_calcs.append(self.int_dict)

		loadModString = 'AFLOWpi.scfuj.transport_prep(oneCalc,ID)'
		self.addToAll(block='PREPROCESSING',addition=loadModString)		
		self.change_input('&control','calculation','"scf"')

		self.tight_banding=False
		


		for temp_index in range(len(temperature)):
			# two flags for whether or not to run the scf/nscf/pdos
			# if doing the calculation at more than one temperature.
			# transport calcs for each temp are considered separate
			# steps in the workflow.
			run_trans_prep=True
			if self.tight_banding==True:
				run_trans_prep=False
			


			run_scf=True
			if self.scf_complete==True:
				run_scf=False



			self.tight_banding=True
			self.scf_complete=True
			

			#only spawn off a new step if it's not the first temperature in the list
			if temp_index!=0:
				self.new_step(update_positions=True,update_structure=True)
				self.initial_calcs.append(self.int_dict)
			loadModString = "oneCalc,ID = AFLOWpi.scfuj.run_transport(__submitNodeName__,oneCalc,ID,temperature=%s,run_scf=%s,run_transport_prep=%s,epsilon=%s,run_bands=%s)"%(temperature[temp_index],run_scf,run_trans_prep,epsilon,run_bands)
			self.addToAll(block='RUN',addition=loadModString)		

			calc_type='Optical and Transport Properties at %sK' % temperature[temp_index]
			print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)

			## no temperature parameter for WanT bands so only run 
			## it once if run_bands=True in the input the method.
			run_bands=False



	def scfuj(self,thresh=0.1,nIters=20, paodir=None,relax='scf',mixing=0.0):
		'''
		Wrapper method to call AFLOWpi.scfuj.scfPrep and AFLOWpi.scfuj.run in the high level 
		user interface. Adds a new step to the workflow.
		

		Arguments:
		      self: the _calcs_container object

		Keyword Arguments:
		      thresh (float): threshold for self consistent hubbard U convergence
		      niters (int): max number of iterations of the acbn0 cycle
		      paodir (string): the path of the PAO directory. This will override 
		                        an entry of the paodir in the AFLOWpi config file 
		                        used for the session
		      mixing (float): the amount of the previous acbn0 U iteration to mix into
                                      the current (only needed when there is U val oscillation)
		Returns:
		      None

		'''		

		self.scf_complete=True
		self.tight_banding=True
		self.load_index+=1
		self.type='scfuj'
		self.new_step(update_positions=True,update_structure=True)
		self.initial_calcs.append(self.int_dict)
		self.change_input('&control','calculation','"%s"'%relax)


		self.int_dict = AFLOWpi.scfuj.scfprep(self.int_dict,paodir=paodir)
		AFLOWpi.scfuj.run(self.int_dict,uThresh=thresh, nIters=nIters,mixing=mixing)
		self.initial_calcs.append(self.int_dict)

		calc_type='ACBN0 Self-Consistent Hubbard U'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)
		
	
	def crawl_min(self,mult_jobs=False,grid_density=10,initial_variance=0.02,thresh=0.01,constraint=None,final_minimization='relax'):
		'''
		Wrapper method to call AFLOWpi.pseudo.crawlingMinimization in the high level user interface.
		Adds a new step to the workflow.

		Arguments:
		      self: the _calcs_container object

		Keyword Arguments:
		      mult_jobs (bool): if True split the individual scf jobs into separate cluster jobs 
		                     if False run them serially
		      grid_density (int): controls the number of calculations to generate for the minimization
		                           num_{calcs}=grid_density^{num_{parameters}-num_{constraints}}
		      initial_variance (float): amount to vary the values of the parameters from the initial
		                                 value. i.e. (0.02 = +/-2% variance)
		      thresh (float): threshold for $\DeltaX$ of the lattice parameters between brute force
		                  minimization iterations.
           	      constraint (list): a list or tuple containing two entry long list or tuples with
  	                                  the first being the constraint type and the second the free 
				          parameter in params that its constraining for example in a 
				          orthorhombic cell: constraint=(["volume",'c'],) allows for A and B
				          to move freely but C is such that it keeps the cell volume the same
				          in all calculations generated by the input oneCalc calculation.		 
		      final_minimization (str): calculation to be run at the end of the brute force minimization
		                                 options include "scf", "relax", and "vcrelax"
				     

		Returns:
		      None

		'''

		self.scf_complete=True
		self.tight_banding==False
		self.load_index+=1
		self.type='crawl_min'
		self.new_step(update_positions=True,update_structure=True)

		self.int_dict=AFLOWpi.pseudo.crawlingMinimization(self.int_dict,mult_jobs=mult_jobs,grid_density=grid_density,initial_variance=initial_variance,thresh=thresh,constraint=constraint)
		self.initial_calcs.append(self.int_dict)

		calc_type='Crawling brute force cell optimization'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)


        def evCurve_min(self,pThresh=25,final_minimization='relax'):
		self.scf_complete=True
		self.type='evCurve_min'
                self.new_step(update_positions=True,update_structure=True)
		self.load_index+=1
		self.int_dict=AFLOWpi.scfuj.evCurveMinimize(self.int_dict,pThresh=pThresh,final_minimization=final_minimization)
		self.initial_calcs.append(self.int_dict)


        def dos(self,kpFactor=2,project=True):
		'''
		Wrapper method to call AFLOWpi.prep.doss in the high level user interface.
		Adds a new step to the workflow.

		Arguments:
		      self: the _calcs_container object

		Keyword Arguments:
		      kpFactor (float): factor to which the k-point grid is made denser in each direction
		      project (bool): if True: do the projected DOS after completing the DOS

		Returns:
		      None

		'''
		self.tight_banding=False
		self.type='dos'
		self.new_step(update_positions=True,update_structure=True)
		self.int_dict = doss(self.int_dict,kpFactor=kpFactor)
		postfix = ConfigSectionMap('run','execpostfix')

		if len(re.findall(r'npool',postfix))!=0 or len(re.findall(r'nk',postfix))!=0:
			self.change_input('&control','wf_collect','.TRUE.')#,change_initial=False)		

		AFLOWpi.run.scf(self.int_dict)
		AFLOWpi.run.dos(self.int_dict)	

		if project==True:
			AFLOWpi.run.pdos(self.int_dict)

		calc_type='Density of States'
		if project==True:
			calc_type+=' with atomic projection'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)


        def bands(self,dk=None,nk=100):
		'''
		Wrapper method to write call AFLOWpi.prep.bands for calculating the Electronic Band
		Structure. 

		Arguments:
		      calcs (dict): a dictionary of dicionaries representing the set of calculations

		Keyword Arguments:
		      dk (float): the density in the Brillouin zone of the k point sampling along the
	                           entirety of the path between high symmetry points.
		      nk (int): the approximate number of sampling points  in the Brillouin Zone along
 	                         the entirety of the path between high symmetry points. Points are 
			         chosen so that they are equidistant along the entirety of the path.
			         The actual number of points will be slightly different than the 
			         inputted value of nk. nk!=None will override any value for dk.

		Returns:
		      None

		'''
		self.tight_banding=False
		self.type='bands'
		self.new_step(update_positions=True,update_structure=True)

		postfix = ConfigSectionMap('run','execpostfix')
		if len(re.findall(r'npool',postfix))!=0 or len(re.findall(r'nk',postfix))!=0:
			self.change_input('&control','wf_collect','.TRUE.')#,change_initial=False)		
	
		AFLOWpi.run.scf(self)

		self.int_dict = bands(self.int_dict,dk=dk,nk=nk)
		AFLOWpi.run.bands(self)

		calc_type='Electronic Band Structure'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)


	def pseudo_test_brute(self,ecutwfc,dual=[],sampling=[],conv_thresh=0.01,constraint=None,initial_relax=None,
			      min_thresh=0.01,initial_variance=0.05,grid_density=7,mult_jobs=False,options=None):


		self.load_index+=1
		self.type='brute_pseudotest'
		self.new_step(update_positions=True,update_structure=True)
		self.initial_calcs.append(self.int_dict)		
#		self.initial_calcs.append(self.int_dict)

		AFLOWpi.pseudo.brute_test(self.int_dict,ecutwfc,dual=dual,sampling=sampling,constraint=None,thresh=conv_thresh,initial_variance=initial_variance,grid_density=grid_density,mult_jobs=mult_jobs,)
		calc_type='Brute Force Pseudotesting'
		print '\nADDING STEP #%02d: %s'% (self.step_index,calc_type)

	

#		if test_type=="phonon":
#			pass
#		if test_type=="all_electron":
#			pass


	def __addToInit(self,block=None,addition=None):
		if block==None or addition==None:
			return
		for ID,oneCalc in self.initial_calcs[0].iteritems():
			AFLOWpi.prep.addToBlockWrapper(oneCalc,ID,block,addition)    

	def submit(self):
		if self.submit_flag==False:
			try:
				init_calcs=initial_calcs=self.initial_calcs[0]
				workflow_add_command='''oneCalc['_AFLOWPI_WORKFLOW_']=%s'''%self.workflow
				

			except:
				print 'Must add at least one step to workflow. Exiting'
				raise SystemExit
			AFLOWpi.run.addatexit__(self.__addToInit,block='PREPROCESSING',addition=workflow_add_command)
			AFLOWpi.run.addatexit__(AFLOWpi.run.submitFirstCalcs__,init_calcs,)
		self.submit_flag=True

class plotter:
	'''
	Class for adding common plotting functions from AFLOWpi.plot module to the high level user 
	interface. 

	'''
	def __init__(self,calcs):
		self.calcs=calcs
	
	# def epsilon(self,nm=False,runlocal=False):
	# 	'''
	# 	Wrapper method to call AFLOWpi.plot.epsilon in the high level user interface.

	# 	Arguments:
	# 	      - self -- the plotter object

	# 	Keyword Arguments:
	# 	      - nm -- whether to plot in nanometers for spectrum or eV for energy
	# 	      - runlocal -- a flag to choose whether or not to run the wrapped function now
	#                 or write it to the _ID.py to run during the workflow
		
	# 	Returns:
	# 	      None

	# 	'''

	# 	AFLOWpi.plot.epsilon(self.calcs,nm=nm,runlocal=runlocal)

	# 	calc_type='Plot Absorbtion'
	# 	print '                 %s'% (calc_type)



	def transport(self,runlocal=False,postfix=''):
		'''
		Wrapper method to call AFLOWpi.plot.epsilon in the high level user interface.

		Arguments:
		      self: the plotter object

		Keyword Arguments:
		      nm (bool): whether to plot in nanometers for spectrum or eV for energy
		      runlocal (bool): a flag to choose whether or not to run the wrapped function now
	                                or write it to the _ID.py to run during the workflow
		
		Returns:
		      None

		'''

		AFLOWpi.plot.transport_plots(self.calcs,runlocal=runlocal,postfix=postfix)
#		AFLOWpi.plot.kappa(self.calcs,runlocal=runlocal)
#		AFLOWpi.plot.seebeck(self.calcs,runlocal=runlocal)
#		AFLOWpi.plot.sigma_seebeck(self.calcs,runlocal=runlocal)
#		AFLOWpi.plot.ZetaT(self.calcs,runlocal=runlocal)
#		AFLOWpi.plot.conduction(self.calcs,runlocal=runlocal)
		

		calc_type='Plot Optical and Transport properties'
		print '                 %s'% (calc_type)

	def opdos(self,yLim=[-10,10],LSDA=False,runlocal=False,postfix=''):
		'''
		Wrapper method to call AFLOWpi.plot.opdos in the high level user interface.

		Arguments:
		      self: the plotter object

		Keyword Arguments:
		      yLim (list): a tuple or list of the range of energy around the fermi/Highest
		                    occupied level energy that is to be included in the plot.
		      LSDA (bool): Plot the up and down of a spin polarized orbital projected DOS
		                    calculation.
		      runlocal (bool): a flag to choose whether or not to run the wrapped function now
               	                        or write it to the _ID.py to run during the workflow
		      postfix (string): a string of an optional postfix to the plot filename for every
		                         calculation.

		Returns:
		      None

		'''

		AFLOWpi.plot.opdos(self.calcs,yLim=yLim,LSDA=LSDA,runlocal=runlocal,postfix=postfix)

		calc_type='Plot Orbital Projected DOS'
		print '                 %s'% (calc_type)


	def phonon(self,runlocal=False,postfix='',THz=True):
		AFLOWpi.plot.phonon(self.calcs,runlocal=runlocal,postfix=postfix,THz=THz)

		calc_type='Plot Phonon Bands and DOS'

		print '                 %s'% (calc_type)		

	def bands(self,yLim=[-10,10],DOSPlot='',LSDA=False,runlocal=False,postfix=''):
		'''
		Wrapper method to call AFLOWpi.plot.bands in the high level user interface.

		Arguments:
		      self: the plotter object

		Keyword Arguments:
		      yLim (list): a tuple or list of the range of energy around the fermi/Highest
		                occupied level energy that is to be included in the plot.
		      DOSPlot (str): a string that flags for the option to have either a DOS plot
		                   share the Y-axis of the band structure plot. 

				   Options include:
				   ""      | A blank string (default) will cause No Density of
				           | States plotted alongside the Band Structure

				   "APDOS" | Atom Projected Density of States
				   "DOS"   | Normal Density of States

				   
		      
		      LSDA (bool): Plot the up and down of a spin polarized orbital projected DOS
		                calculation.
		      runlocal (bool): a flag to choose whether or not to run the wrapped function now
               	                    or write it to the _ID.py to run during the workflow
		      postfix (str): a string of an optional postfix to the plot filename for every
		                   calculation.

		Returns:
		      None

		'''

		AFLOWpi.plot.bands(self.calcs,yLim=yLim,DOSPlot=DOSPlot,LSDA=LSDA,runlocal=runlocal,postfix=postfix)

		calc_type='Plot Electronic Band Structure'
		if DOSPlot=='DOS':
			calc_type+=' with Density of States'
		if DOSPlot=='APDOS':
			calc_type+=' with Atom Projected Density of States'
		print '                 %s'% (calc_type)

	def dos(self,yLim=[-10,10],LSDA=False,runlocal=False,postfix=''):
		'''
		Wrapper method to call AFLOWpi.plot.dos in the high level user interface.

		Arguments:
		      self: the plotter object

		Keyword Arguments:
		      yLim (list): a tuple or list of the range of energy around the fermi/Highest
		                    occupied level energy that is to be included in the plot.
		      LSDA (bool): Plot the up and down of a spin polarized DOS
		                    calculation.
		      runlocal (bool): a flag to choose whether or not to run the wrapped function now
               	                 or write it to the _ID.py to run during the workflow
		      postfix (str): a string of an optional postfix to the plot filename for every
		                      calculation.

		Returns:
		      None

		'''

		AFLOWpi.plot.dos(self.calcs,yLim=yLim,LSDA=LSDA,runlocal=runlocal,postfix=postfix)

		calc_type ='Plot Density of States'
		print '                 %s'% (calc_type)





#########################################################################################################################
#########################################################################################################################
##MODIFY QE INPUT
#########################################################################################################################
######################################################################################################################### 

def __num_bands(oneCalc):
    '''
    attempt to find the number of kohn-sham orbitals after the scf calculation to find
    the value of the "nbnd" parameter in the "&system" namelist of QE input files.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation

    Returns:
          The value for the nbnd parameter in QE nscf calculations
    
    '''

    subdir = oneCalc['_AFLOWPI_FOLDER_']

    for old_ID in oneCalc['prev']:
	try:
		oldFile = os.path.join(subdir,'%s.out' % old_ID)
		if os.path.exists(oldFile):
			with open(oldFile,'r') as outfileObj:
				outfile=outfileObj.read()
				match1 = re.search(r'Kohn-Sham states=',outfile)
				m1 = match1.end()
				m2 = m1 + 15
				nbnd = int(2.0*int(outfile[m1:m2]))
				print 'Number of bands to be Calculated %s: '% nbnd
				logging.info('Number of bands to be Calculated %s: '% nbnd)
				break
	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)

    return nbnd


def doss(calcs,kpFactor=1.5):
	'''
	Wrapper function to write the functio n AFLOWpi.prep.__oneDoss to the _ID.py

	Arguments:
	     calcs (dict): a dictionary of dicionaries representing the set of calculations 

	Keyword Arguments:
	     kpFactor (float): the factor to which we make each direction in the kpoint grid denser
	      
	Returns:
	      The identical "calcs" input variable

	'''

	loadModString = 'AFLOWpi.prep.__oneDoss(oneCalc,ID,kpFactor=%s)'%kpFactor
	__addToAll(calcs,block='PREPROCESSING',addition=loadModString)
	'''get the fermi level from here'''
	AFLOWpi.prep.__addToAll(calcs,block='POSTPROCESSING',addition='AFLOWpi.retr.__writeEfermi(oneCalc,ID)')
	return calcs


import math
@newstepWrapper(__check_lock)
def __oneDoss(oneCalc,ID,kpFactor=1.5):
    '''
    Converts an scf calculation to an nscf for calculating DOS.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          kpFactor (float): the factor to which we make each direction in the kpoint grid denser
    
    Returns:
          the one calculation dictionary object with oneCalc['_AFLOWPI_INPUT_'] converted to a 
	  DOS nscf calcuation.
    
    '''

    try:
        logging.debug('entering __oneDoss')
        d=oneCalc

        f=ID
	output_calcs = {}
        
        subdir = d['_AFLOWPI_FOLDER_']
        inputfile = d['_AFLOWPI_INPUT_']
        

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)


    inputDict=AFLOWpi.retr.__splitInput(inputfile)
    inputDict['&control']['calculation']="'nscf'"

    nbnd = AFLOWpi.prep.__num_bands(oneCalc)

    try:
	    '''adds nbnd to input file'''
	    inputDict['&system']['nbnd']=str(nbnd)

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    try:
        scfFileString=oneCalc['_AFLOWPI_INPUT_']
	mod = inputDict['K_POINTS']['__modifier__'].upper()
	mod=mod.strip('{}()')
	if mod=='GAMMA':
		inputDict['K_POINTS']['__content__']='2 2 2 0 0 0'
		inputDict['K_POINTS']['__modifier__']='{automatic}'
	else:
		scfKPointString =  AFLOWpi.qe.regex.k_points(scfFileString)
		scfKPointSplit = scfKPointString.split()
		scfKPointSplit = [float(x) for x in scfKPointString.split()]
		for kpoint in range(len(scfKPointSplit)-3):
		    scfKPointSplit[kpoint] = str(int(math.ceil(scfKPointSplit[kpoint]*kpFactor)))+' '
		for kpoint in range(3,len(scfKPointSplit)):
		    scfKPointSplit[kpoint] = '0 '
		newKPointString = ''.join(scfKPointSplit)

		inputDict['K_POINTS']['__content__']=newKPointString
		inputDict['K_POINTS']['__modifier__']='{automatic}'

        inputfile = AFLOWpi.retr.__joinInput(inputDict)

    except Exception,e:                    
        AFLOWpi.run._fancy_error_log(e)

    try:
	    
        a = ID+'.in'
        new_inputfile = open(os.path.join(subdir,a),'w')
        new_inputfile.write(inputfile)
        new_inputfile.close()

        d['execFile']= os.path.join(subdir,'_%s.py' % ID)

        output_calcs = d
        output_calcs['_AFLOWPI_INPUT_'] = inputfile
        '''set prev to current ID so we can check if we've already done the transform'''
        try:
            output_calcs['prev'].append(f)
        except:
            output_calcs['prev']=[f]
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)


    logging.debug('exiting __oneDoss')                                              
    return output_calcs,ID


def modifyInputPrefixPW(calcs,pre):
    '''
    A Wrapper function that is used to write a function to the __ID.py

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations 

    Returns:
          None

    '''
	
    AFLOWpi.prep.__addToAll(calcs,block='PREPROCESSING',addition="oneCalc,ID = AFLOWpi.prep.__modifyInputPrefixPW(oneCalc,ID)")    
    return calcs

import __main__

def __modifyInputPrefixPW(oneCalc,ID):
    '''
    Assigns the new prefix to a the calculation inputs in the form _ID
	
    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
	  ID (str): The ID string for the particular calculation and the value of the the
                  prefix without the leading "_" in front.

    Returns:
          None
	  
    '''

    try:

	    inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
	    new_prefix='_'+ID
	    inputDict['&control']['prefix']=repr(new_prefix)
	    inputString = AFLOWpi.retr.__joinInput(inputDict)
	    oneCalc['_AFLOWPI_INPUT_']=inputString
	    oneCalc['_AFLOWPI_PREFIX_']=new_prefix
	    AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
	    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as newIn:
		    newIn.write(inputString)
    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)


    return oneCalc,ID
			
			
			
    
def modifyNamelistPW(calcs,namelist,parameter,value,runlocal=False):
    '''
    A Wrapper function that is used to write the function AFLOWpi.prep.__modifyNameListPW 
    to the _ID.py. If the value is intended to be a string in the QE input file, it must
    be dually quoted i.e. value="'scf'" will become 'scf' in the input file.

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations 
          namelist (str): a string of the fortran namelist that the parameter is in
          parameter (str): a string of the parameter name
          value: the value of that parameter 

    Keyword Arguments:
          runlocal (bool): a flag to choose whether or not to run the wrapped function now
	                or write it to the _ID.py to run during the workflow.

    Returns:
          Either the identical set of calculations if runlocal == False or the set of
          calculations with the parameter's value changed in their oneCalc['_AFLOWPI_INPUT_']
          if runlocal==True

    '''
 
    if runlocal==True:
	    for ID,oneCalc in calcs.iteritems():
		    temp,temp_ID=__modifyNamelistPW(oneCalc,ID,namelist,parameter,value)
		    calcs[ID]=temp
    else:
	    addit="oneCalc,ID = AFLOWpi.prep.__modifyNamelistPW(oneCalc,ID,'%s','%s',%s)" %(namelist,parameter,repr(value))
	    AFLOWpi.prep.__addToAll(calcs,block='PREPROCESSING',addition=addit)
    return calcs


def __modifyNamelistPW(oneCalc,ID,namelist,parameter,value):
    '''
    A Wrapper function that is used to write the function AFLOWpi.prep.__modifyNameListPW 
    to the _ID.py. If the value is intended to be a string in the QE input file, it must
    be dually quoted i.e. value="'scf'" will become 'scf' in the input file. If the value
    is equal to None (without quotes) it will remove that parameter from the namelist.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step
          namelist (str): a string of the fortran namelist that the parameter is in
          parameter (str): a string of the parameter name
          value: the value of that parameter 

    Returns:
          oneCalc object representing the calculation but with the paramater's value
          changed or added in oneCalc['_AFLOWPI_INPUT_'] and the same ID as the input to 
	  the function.
          
    '''
    
    inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
    
    try:
        if value==None:
		try:
			del inputDict[namelist][parameter]
	        except:
			pass
	else:
		inputDict[namelist][parameter]=value
	 
        inputString = AFLOWpi.retr.__joinInput(inputDict)
        oneCalc['_AFLOWPI_INPUT_']=inputString
        AFLOWpi.prep.__saveOneCalc(oneCalc,ID)
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as newIn:
            newIn.write(inputString)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        logging.warning('namelist %s and/or parameter %s do not exist in input file.' % (namelist,parameter))

    return oneCalc,ID



def lockAtomMovement(calcs):
    '''
    A Wrapper function that writes the function AFLOWpi.prep.__freezeAtoms to the __ID.py

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations 

    Returns:
          None

    '''

    for ID,oneCalc in calcs.iteritems():
        __addToBlock(oneCalc,ID,'PREPROCESSING','''oneCalc = AFLOWpi.prep.__freezeAtoms(oneCalc,ID) ''')

def unlockAtomMovement(calcs):
    '''
    A Wrapper function that writes the function AFLOWpi.prep.__unfreezeAtoms to the __ID.py

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations 

    Returns:
          None

    '''

    for ID,oneCalc in calcs.iteritems():
        __addToBlock(oneCalc,ID,'PREPROCESSING','''oneCalc = AFLOWpi.prep.__unfreezeAtoms(oneCalc,ID)''')
            
            

def __freezeAtoms(oneCalc,ID):
    '''
    Modifies a QE input file to have all the atom movement flags in the ATOMIC_POSITIONS 
    card be set to 0 which means they cannot move during a ionic or variable cell relax
    calculation.
    
    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step


    Returns:
          oneCalc object representing the calculation but with the flags for atom
	  movement all set to 0.

    '''

    inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
    atomPos = inputDict['ATOMIC_POSITIONS']['__content__']
    positions,flags = AFLOWpi.retr.detachPosFlags(atomPos)
    flagArray=[]
    for item in range(len(positions.split('\n'))):
        flagArray.append(' 0 0 0')

    frozenFlags='\n'.join(flagArray)

    newPosString = AFLOWpi.retr.attachPosFlags(positions,frozenFlags)
    inputDict['ATOMIC_POSITIONS']['__content__']=newPosString

    newInput = AFLOWpi.retr.__joinInput(inputDict)
    oneCalc['_AFLOWPI_INPUT_']=newInput
    AFLOWpi.prep.__saveOneCalc(oneCalc,ID)

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as inputFile:
        inputFile.write(newInput)

    return oneCalc

def __unfreezeAtoms(oneCalc,ID):
    '''
    Modifies a QE input file to have all the atom movement flags in the ATOMIC_POSITIONS 
    card be set to 1 which means they are allowed to move during a ionic or variable cell
    relax calculation which is default in QE.
    
    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step


    Returns:
          oneCalc (dict): object representing the calculation but with the flags for atom
	  movement all set to 1.

    '''

    inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
    atomPos = inputDict['ATOMIC_POSITIONS']['__content__']
    positions,flags = AFLOWpi.retr.detachPosFlags(atomPos)
    flagArray=[]
    for item in range(len(positions.split('\n'))):
        flagArray.append(' 1 1 1')

    frozenFlags='\n'.join(flagArray)

    newPosString = AFLOWpi.retr.attachPosFlags(positions,frozenFlags)
    inputDict['ATOMIC_POSITIONS']['__content__']=newPosString

    newInput = AFLOWpi.retr.__joinInput(inputDict)
    oneCalc['_AFLOWPI_INPUT_']=newInput
    AFLOWpi.prep.__saveOneCalc(oneCalc,ID)

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as inputFile:
        inputFile.write(newInput)

    return oneCalc




def changeCalcs(calcs,keyword='calculation',value='scf'):
    '''
    A Wrapper function that writes the function AFLOWpi.prep.__changeCalcs to the __ID.py

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations

    Keyword Arguments:
	  keyword (str): a string which signifies the type of change that is to be made
	  value: the value of the choice.

    Returns:
          The identical set of calculations as the input to this function

    '''

    loadModString = "oneCalc,ID = AFLOWpi.prep.__oneChangeCalcs(oneCalc,ID,keyword='%s',value=%s)" %(keyword,repr(value))
    AFLOWpi.prep.__addToAll(calcs,block='PREPROCESSING',addition=loadModString)
    return calcs


def updateStructs(calcs,update_structure=True,update_positions=True):
    '''
    A Wrapper function that writes the function AFLOWpi.prep.__oneUpdateStructs to the __ID.py

    Arguments:
          calcs (dict): a dictionary of dicionaries representing the set of calculations

    Keyword Arguments:
	  update_structure (bool): if True update the cell parameter if possible from the
	                        output of previous calculations in the workflow.
	  update_positions (bool): if True update the atomic positions if possible from the
	                        output of previous calculations in the workflow.

    Returns:
          The identical set of calculations as the input to this function

    '''

    loadModString ="oneCalc,ID = AFLOWpi.prep.__oneUpdateStructs(oneCalc,ID,update_structure=%s,update_positions=%s)" %(update_structure,update_positions)

    AFLOWpi.prep.__addToAll(calcs,block='PREPROCESSING',addition=loadModString)
    return calcs




import numpy as np    
@newstepWrapper(__check_lock)
def __oneUpdateStructs(oneCalc,ID,update_structure=True,update_positions=True,override_lock=False):
    '''
    Attempts to read output files from previous steps in the workflow and update the input file of
    the current step with updated cell parameters and/or atomic positions that were calculated in
    previous steps. 

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
	  update_structure (bool): if True update the cell parameter if possible from the
	                        output of previous calculations in the workflow.
	  update_positions (bool): if True update the atomic positions if possible from the
	                        output of previous calculations in the workflow.

	  
	  override_lock (bool): <DEFUNCT OPTION: CONSIDER FOR REMOVAL>
	                         override the __check_lock function wrapped around __oneUpdateStructs

    Returns:
          calculations with their cell parameters and/or atomic positions updated from the last 
	  variable cell or ionic relaxation calculation. If there is none previously then the
	  parameters and position remain the same.
    
    
    '''

    try:
        logging.debug('entering __oneUpdateStructs')
        d=oneCalc
        f=ID
        output_calcs = {}
        
        subdir = d['_AFLOWPI_FOLDER_']
        inputfile = d['_AFLOWPI_INPUT_']
        inputDict=AFLOWpi.retr.__splitInput(inputfile)
        for prevOut in reversed(oneCalc['prev']):
            oldFile = os.path.join(subdir,'%s.out' % prevOut)
            if not os.path.exists(oldFile):
                continue
            try:
                outputfile = file(oldFile,'r').read()
                cellParaRE=AFLOWpi.qe.regex.cell_parameters('','content','regex')
                cellPara = cellParaRE.findall(outputfile)

                if len(cellPara) != 0:
                    cellPara=cellPara[-1]
                    #Get alat
                    alatRE = re.compile('alat=\s*(\d+.\d+)')
                    alatSearch = re.compile(r'(?:CELL_PARAMETERS)\s*.*alat[\D]*([0-9.]*)',re.M)
                    alatFind = alatSearch.findall(outputfile)

                    alat = float(alatFind[-1])

                    #Get ibrav
                    splitInput = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
                    try:
                        ibrav=int(splitInput['&system']['ibrav'])
                    except:
                        ibrav=0
                    #Make cell parameter matrix
                    temp = []
                    splitPara = cellPara.split('\n')
                    for item in splitPara:
                        if len(item)!=0:
                            temp.append([float(x) for x in item.split(' ') if len(x)!=0])

                    cellParaMatrix = alat*np.array(temp).astype(np.float)
                    #Update cell dimensions according to ibrav and alat
                    ibravDict = AFLOWpi.retr.__free2celldm(cellParaMatrix,ibrav=ibrav)                    
                    ibravDict['ibrav']=ibrav

                    break
                else:
                    pass

            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)

        for prevOut in reversed(oneCalc['prev']):
            oldFile = os.path.join(subdir,'%s.out' % prevOut)
            if not os.path.exists(oldFile):
                continue
            outputfile = file(oldFile,'r').read()
            try:
                atmPos  = AFLOWpi.qe.regex.atomic_positions(outputfile)


                if len(atmPos.strip())!=0:
                    break
            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)


        splitInput = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])   
        try:
            if update_structure==True:

                for item in ibravDict.items():
                    splitInput['&system'].update({item[0]:item[1]})
        except:
            pass
        try:
            if update_positions==True:
                if len(atmPos.strip())!=0:
                    atom_pos_input = splitInput['ATOMIC_POSITIONS']['__content__']
                    fakeFlags,flags=AFLOWpi.retr.detachPosFlags(atom_pos_input)
                    atmPos,outputflags=AFLOWpi.retr.detachPosFlags(atmPos)
                    atmPos=AFLOWpi.retr.attachPosFlags(atmPos,flags)
                    splitInput['ATOMIC_POSITIONS']['__content__']=atmPos        
                    splitInput['ATOMIC_POSITIONS']['__modifier__']='{crystal}'
        except Exception,e:
            pass
#            AFLOWpi.run._fancy_error_log(e)

        inputfile = AFLOWpi.retr.__joinInput(splitInput)

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    try:
        '''Update CELL_PARAMETERS and ATOMIC_POSITIONS from vc-relax or relax output'''
        output_calcs = copy.deepcopy(d)
        output_calcs['_AFLOWPI_INPUT_'] = inputfile
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)


    logging.debug('exiting __oneUpdateStructs')
    '''written funny because output needs to be oneCalc,ID'''
    try:
	    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as newIn:
		    newIn.write(output_calcs['_AFLOWPI_INPUT_'])

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    return output_calcs,ID


@newstepWrapper(__check_lock)
def __oneChangeCalcs(oneCalc,ID,keyword='calculation',value='scf'):
    '''
    ############
    ##OBSOLETE##
    ############
    Changes a PWSCF calculation's input depending on the option chosen for input variable 
    "keyword". 

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          keyword (str): the type of change that is to be made to the input file. 
	              Options: "calculation","ecutwfc","pseudos","kpoints"

	  value: the value of the choice.
	             Options with possible values:

      	             "calculation" | 'scf','vc-relax' or 'relax'
		     "ecutwfc"     | integer value
		     "pseudos"     | a string of a local directory containing PP files
		     "kpoints"     | a string containing some kind of list of k points

		      
		    
    Returns:
          returns oneCalc with oneCalc['_AFLOWPI_INPUT_'] modified in some way depending on 
	  the option chosen for "keyword" and the ID label for calculation at that step
	  in the workflow.

    '''

    try:
        """
        have to keep the index for tracking purposes between _<ID>.py scripts
        """
        logging.debug('entering __oneChangeCalcs')
        d=oneCalc
	prefix = oneCalc['_AFLOWPI_PREFIX_'][1:]
        f=ID



        output_calcs = {}
        subdir = d['_AFLOWPI_FOLDER_']
        oldFile = os.path.join(subdir,'%s.out' % oneCalc['prev'][-1])
 
        if not os.path.exists(oldFile):
            logging.debug('TEMP FILE DOES NOT EXIST FOR __oneChangeCalcs...EXITING')
#            raise SystemExit


    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    try:

	logging.info("Updating %s of %s to %s"%(keyword,prefix,value))
	print "Updating %s of %s to %s"%(keyword,prefix,value)

	if keyword == 'calculation':
	        '''we need to change calculation to vc-relax'''
	        inputfile = d['_AFLOWPI_INPUT_']

#	        replaceCalcRE = re.compile(r"calculation\s*=\s*'\s*(.+?)\s*'")
#		replaceCalc = replaceCalcRE.findall(inputfile)[0]
		# if replaceCalc == 'vc-relax' or replaceCalc == 'relax':
		# #try to uppdate Crystal structure
                #     logging.info("Updating Crystal Structure of %s from %s calculation"%(prefix,replaceCalc))
                #     print "Updating Crystal Structure of %s from %s calculation"%(prefix,replaceCalc)

#                d = __oneUpdateStructs(oneCalc, oneCalc['prev'],update_structure=update_structure,update_positions=update_positions)

         		#Change calculation type
                inputfile = d['_AFLOWPI_INPUT_']
                inputDict=AFLOWpi.retr.__splitInput(inputfile)                    
                inputDict['&control']['calculation']=value
                inputfile=AFLOWpi.retr.__joinInput(inputDict)


	elif keyword == 'kpoints':
		replaceKpts = re.compile(r"K_POINTS.*?\n(\d+\s\d+\s\d+)")
		inputfile = d['_AFLOWPI_INPUT_']
                #Change kpoint grid
                inputfile = replaceKpts.sub("K_POINTS {automatic}\n%s"%value,inputfile)

	elif keyword == 'ecutwfc':
		replaceRegEx = re.compile(r"ecutwfc\s*?=.*?\n")
                inputfile = d['_AFLOWPI_INPUT_']
                #Change energy cutoff
                inputfile = replaceRegEx.sub("ecutwfc = %s,\n"%value,inputfile)
		if re.search(r"ecutrho",inputfile):
			replaceRegEx = re.compile(r"ecutrho\s*?=.*?\n")
			#Change density cutoff
			inputfile = replaceRegEx.sub("ecutrho = %f,\n"%(float(value)*4.0),inputfile)

	elif keyword == 'pseudos':
		inputfile = d['_AFLOWPI_INPUT_']
		new_pseudodir = value
		for key in d:
        	        v = d[key]
                        if re.search(r'_AFLOWPI_[A-Z][0-9]*_', key):
                                speciesRe = ''.join([i for i in v if not i.isdigit()])

                                vp = __getPseudofilename(speciesRe,new_pseudodir)

                                try:
                                        a = os.path.join(new_pseudodir,vp)
                                	b = os.path.join(subdir, vp)
                                	if not os.path.exists(b):
	                                    shutil.copy(a, b)
                                except AttributeError:
                                        logging.debug('Cannot find pseudopotential files in %s ...Exiting' % new_pseudodir)
                                        print 'cannot find correct PAO files in %s ...Exiting' % new_pseudodir
                                        raise SystemExit
         
                 
                                #Replace PP name in inputfile
                 
                                ppLineRE=re.compile("%s.*\s+\S*\.UPF"%v) # v is element symbol
                                ppLine = ppLineRE.findall(inputfile)
                                replacePpRE = re.compile("\S*\.UPF")
                                replacePp = replacePpRE.sub(vp,ppLine[0])

                                inputfile = ppLineRE.sub(replacePp,inputfile)

	else:
		replaceRegEx = re.compile(r"%s\s*?=.*?\n"%keyword)
		inputfile = d['_AFLOWPI_INPUT_']
                #Change value of keyword
                inputfile = replaceRegEx.sub("%s = %s,\n"%value,inputfile)

    except Exception,e:
        print e
	AFLOWpi.run._fancy_error_log(e)

    try:

        calc_label = ID
	a = calc_label+'.in'
	
        new_inputfile = open(os.path.join(subdir,a),'w')
        new_inputfile.write(inputfile)
        new_inputfile.close()

        d['execFile']= os.path.join(subdir,'_%s.py' % calc_label)

        output_calcs = d
        output_calcs['_AFLOWPI_INPUT_'] = inputfile
        '''set prev to current ID so that we can check if we have already transformed'''
        # try:
        #     output_calcs['prev'].append(f)
        # except:
        #     output_calcs['prev']=[f]
        # """save it so if it's loaded by a scipt it'll have the correct vals"""
        # __saveOneCalc(oneCalc,ID)

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)



    logging.debug('exiting __oneChangeCalcs')
    return output_calcs,calc_label




#########################################################################################################################
#########################################################################################################################
## MISC
######################################################################################################################### 
#########################################################################################################################


def varyCellParams(oneCalc,ID,param=(),amount=0.15,steps=8,constraint=None):
    '''
    Forms and returns a set of calcs with varied cell params must be in A,B,C, 
    and, in degrees,alpha,beta,gamma and then returns it.


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step


    Keyword Arguments:
          param (tuple): the params assoc. with the amount and step. i.e. ('celldm(1)','celldm(3))
          amount (float): percentage amount to be varied up and down. i.e (0.04,0.02,0.01)
          steps  (int): how many steps to within each range. i.e (4,5,7) 
          constraint (list): a list or tuple containing two entry long list or tuples with
  	                  the first being the constraint type and the second the free 
			  parameter in params that its constraining for example in a 
			  orthorhombic cell: constraint=(["volume",'c'],) allows for A and B
			  to move freely but C is such that it keeps the cell volume the same
			  in all calculations generated by the input oneCalc calculation.

    Returns:
          A dictionary of dictionaries representing a new calculation set
			  			  
    '''

    oneCalc=copy.deepcopy(oneCalc)

    constraint_type_list=[]
    constraint_var_list=[]
    if constraint!=None:
        try:
            '''if the length of the first entry in constraints isn't a   '''
            '''[constraint,parameter] combination then put it in an array'''            
            if len(constraint[0])!=2:
                constraint=[list(constraint),]

            for constraints in constraint:
                constraint_type_list.append(constraints[0].lower().strip())
                constraint_var_list.append(constraints[1].lower().strip())
                
        except Exception,e:
            print e
            constraint_type=None
            constraint_var =None
    else:
        constraint_type=None
        constraint_var =None


#    constraint_type='volume'
#    constraint_var='c'
#    rangeArray = build_calcs(tuple(ranges.values()),build_type='product')
    inputList=[]

    inputDict = AFLOWpi.retr.__splitInput(oneCalc['_AFLOWPI_INPUT_'])
    paramDict={}
    param=[x.lower() for x in param]

    for namelist,paramlist in inputDict.iteritems():
        for k,v in inputDict[namelist].iteritems():
            if re.match(r'celldm.*',k.lower()):
                paramDict[k.replace(')','').replace('(','')]=float(v)
            if re.match(r'ibrav.*',k.lower()):
                paramDict[k]=int(v)

    paramDict['cosine']=False
    paramDict['degrees']=True
    abcList = ['a','b','c','alpha','beta','gamma']
    abcDict={}
    abcRes = AFLOWpi.retr.celldm2abc(**paramDict)
    for item in range(len(abcRes)):
        abcDict[abcList[item]]=abcRes[item]

    
    ibrav=paramDict['ibrav']
    if ibrav==1 or ibrav==2 or ibrav==3:
        changeList=['a']
    if ibrav==5:
        changeList=['a','gamma']
    if ibrav==6 or ibrav==7 or ibrav==4:
        changeList=['a','c']
    if ibrav==8 or ibrav==9 or ibrav==10 or ibrav==11:
        changeList=['a','b','c']
    if ibrav==12 or ibrav==13:
        changeList=['a','b','c','alpha']
    if ibrav==14:
        changeList=['a','b','c','alpha','beta','gamma']


    
    try:
        numCalc=reduce(lambda x, y: x*y, steps)    
        if (len(changeList)-len(constraint))<3 and numCalc>1000:
            logging.warning('number of unconstrained variable=%s. number of calcs to be generated=%s. This exceeds the limit set at 1000 in AFLOWpi. If you would like to ignore the limit you can modify numCalc variable in AFLOWpi.prep.varyCellParams.' % (len(changeList)-len(constraint),numCalc))
    except:
        pass
    try:

        newDict=copy.deepcopy(abcDict)
        newDict['ibrav']=ibrav

        vol = AFLOWpi.retr.abcVol(**newDict)

	if ibrav in [4,5]:
		vol*=3
	if ibrav in [3,7,9,11,13]:
		vol*=2
	if ibrav in [2,10]:
		vol*=4

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    productList=[]
    for item in range(len(param)):
        if param[item] in changeList:
            if param[item] not in constraint_var_list:

                productList.append(numpy.linspace(abcDict[param[item]]*(1.0-amount[item]),abcDict[param[item]]*(1.0+amount[item]),steps[item]).tolist())

    modifierList = list(it.product(*productList))
    param=[x for x in param if x in changeList]


    modifierDictList=[]
    
    tempDict=collections.OrderedDict({'ibrav':ibrav})
    


    tempCelldmDict={}   
    for val in modifierList:
        counter=0
        for item in param:


            if item not in constraint_var_list:
                tempDict[item]=val[counter]
                counter+=1



        for item in range(len(constraint_type_list)):
            for some_param in range(len(param)):
                if param[some_param] in constraint_var_list:
                    if constraint_type_list[item]=='length':
                        pass

        for item in range(len(constraint_type_list)):
            for some_param in range(len(param)):
                if param[some_param] in constraint_var_list:
                    if constraint_type_list[item]=='volume':
                        volDiv=1
                        for parameter,value in tempDict.items()[1:]:
                            if parameter not in constraint_var_list:
                                if parameter in ['alpha','beta','gamma']:
                                    volDiv*=numpy.sin(value*(numpy.pi/180.0))
                                if parameter in ['a','b','c']:
                                    volDiv*=value


                        remaining = vol/volDiv
                        
                        if constraint_var_list[item] in ['alpha','beta','gamma']:
                            remaining = numpy.arcsin(remaining)*180.0/numpy.pi

                        tempDict[param[some_param]]=remaining

        tempCelldmTuple=AFLOWpi.retr.abc2celldm(**tempDict)[1:]

        for item in range(len(tempCelldmTuple)):
            tempCelldmDict['celldm(%s)'%(int(item)+1)]=tempCelldmTuple[item]


        for k,v in inputDict['&system'].iteritems():
            if k in tempCelldmDict.keys():
                inputDict['&system'][k]=tempCelldmDict[k]
        inputList.append(AFLOWpi.retr.__joinInput(inputDict))



    aflowkeys = {'project':'GRID_MIN','set':'',}
    intoInit={'PROJECT':'GRID_MIN','SET':'','workdir':oneCalc['_AFLOWPI_FOLDER_'],'config':oneCalc['_AFLOWPI_CONFIG_']}
    AFLOWpi.prep.init(**intoInit)
    calcs = calcFromFile(aflowkeys,inputList,workdir=oneCalc['_AFLOWPI_FOLDER_'])


    newConfigPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'GRID_MIN','AFLOWpi','CONFIG.config')

    config = ConfigParser.RawConfigParser()

    config = ConfigParser.RawConfigParser()
    config.read(newConfigPath)
    config.set('prep', 'workdir', oneCalc['_AFLOWPI_FOLDER_'])
    with open(newConfigPath,'w') as fileWrite:    
        config.write(fileWrite)

    try:
	    for ID_new,oneCalc_new in calcs.iteritems():
		    calcs[ID_new]['__walltime_dict__']=oneCalc['__walltime_dict__']
    except:
	    pass

    return calcs
        
def __incrementFileValue(oneCalc,ID,varName='uValue'):
    '''
    ##########################
    ###NEEDS GENERALIZATION###
    ##########################
    adds +1 to the value of a variable defined in the _ID.py files

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          varName (str): a string of the name of the variable that is to be incremented

    Returns:
          None

    '''

    subdir = oneCalc['_AFLOWPI_FOLDER_']    
    fileName = '_'+ID+'.py'
    
    try:
        with open(os.path.join(subdir,fileName),'r') as inputfile:
            inputFileText = inputfile.read()


#            URegex = re.compile(r'%s = (.+)' % varName)
            URegex = re.compile(r'uValue = (.+)')
            uValue = float(URegex.findall(inputFileText)[-1])
            uValue -= 1
            uValueString = str(uValue)
            URegex2 = re.compile(r'uValue = .+')
            newinputfile = URegex2.sub('uValue = ' + uValueString,inputFileText)
        

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        

    try:
        with open(os.path.join(subdir,'_'+ID+'.py'),'w') as outfileObj:
            outfileObj.write(newinputfile)
        

    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)



def __null__(*args,**kwargs):
    '''
    ################
    ###NOT NEEDED###
    ################
    '''
    return 


def __modifyVarVal(oneCalc,ID,varName='uValue',value=None):
    '''
    Modifies the value of a variable defined in the _ID.py files

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          varName (str): a string of the name of the variable that is to be changed
	  value: the value to change it to in the _ID.py

    Returns:
          None

    '''    

    logging.debug('entering __modifyVarVal')
    if value==None:
        return
    subdir = oneCalc['_AFLOWPI_FOLDER_']    
    fileName = '_'+ID+'.py'
    
    try:
        with open(os.path.join(subdir,fileName),'r') as inputfile:
            inputFileText = inputfile.read()


        URegex = re.compile(r'%s = (.+)' % varName)

        try:
            uValueString = str(value)
        except:
            uValueString = value
            
        URegex2 = re.compile(r'%s\s*=\s*.+' % varName)
        newinputfile = URegex2.sub('%s = %s' % (varName,uValueString),inputFileText)
        
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        

    try:
        with open(os.path.join(subdir,'_'+ID+'.py'),'w') as outfileObj:
            outfileObj.write(newinputfile)
        
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
    logging.debug('exiting __modifyVarVal')


def __write_scratch_meta_data_file(oneCalc,ID):
    '''
    Not used. Delete possibly

    '''	
    temp_dir= AFLOWpi.prep.__get_tempdir()
    file_list=os.walk(temp_dir, topdown=True, onerror=None, followlinks=True)
    info_list=[]
    file_path_list=[]
    for file_info in file_list:

	    for file_name in file_info[2]:
		file_path_list.append(os.path.join(file_info[0],file_name))



    for file_name in file_path_list:
	    mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime = os.stat(file_name)
	    info_list.append([str(x) for x in [file_name,mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime]])

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_local_scratch_meta_data.temp'),'w+') as meta_file:
	    for i in info_list:
		    write_string=' '.join(i)+'\n'
		    meta_file.write(write_string)



def __check_restart(oneCalc,ID):
    '''
    Not used. Delete possibly

    '''
    try:
#        inFileString = AFLOWpi.retr.__getInputFileString(oneCalc,ID)
#        inputDict = AFLOWpi.retr.__splitInput(inFileString)
#        mode = inputDict['&control']['restart_mode']
#        mode = mode.lower()
#        mode = mode.strip()
#        mode = mode.strip('"')
#        mode = mode.strip("'")
	
	if oneCalc['__status__']['Restart']!=0:
            return True
        else:
            return False
    except:
        return False


def askAFLOWpiVars(refAFLOWpiVars):
	"""
	Cycle on the keys of the refAFLOWpiVars dictionary and ask to define them

	Arguments:
	      refAFLOWpiVars (dict): the variables in the ref files that you need to input to run the calculation

	Returns:
	      None

	"""

	for k in refAFLOWpiVars.items():
		if k[1] == None:
			try:
				a = input('define '+k[0]+': ')
				refAFLOWpiVars.update({k[0]:a})	
				if type(k[0]) != tuple: 
					raise ValueError, "Must be a tuple"
			except:
				pass

                        


def __passGlobalVar(varname,value):
    '''
    ############
    ##OBSOLETE##
    ############
    Used to define the value of a global variable in this module's global namesapce 
    usually from another module.
	
    Arguments:
          varname (str): name of the global variable
	  value: value of the global variable
       
    Returns:
          None

    '''

    globals()[varname]=value


import string
def __forceSubmitNodeIP(nodeName):
    '''
    ############
    ##OBSOLETE##
    ############
    Creates the global variable __submitNodeName__ and gives it the value of nodeName


    Arguments:
          nodeName (str): a string containing the name of the submit node

    Returns:
          None
    '''
    
    global __submitNodeName__
    __submitNodeName__ = nodeName


def __announcePrint(string):
    '''
    A debugging tool used to accentuate a string os it can be easily picked out from stdout

    Arguments:
          string (str): a string

    Returns:
          None
    '''

    print '#####################################################################################'
    print string
    print '#####################################################################################'



def generateAnotherCalc(old,new,calcs):
    """
    ############
    ##OBSOLETE##
    ############
    Modify the calculation in each subdir and update the master dictionary

    Arguments:
          old (str): string to replace
          new (str): replacement string
          calcs (dict): dictionary of dictionaries of calculations

    Returns:
          A new set of calculations with a new ID of the hash of the new input strings
    
    """

    new_calcs = copy.deepcopy(calcs)
    output_calcs = {}

    for f,d in new_calcs.iteritems():
	try:
		subdir = new_calcs[f]['_AFLOWPI_FOLDER_']
		a = f+'.in'
		inputfile = open(os.path.join(subdir,a),'r').read()
		inputfile = re.sub(old,new,inputfile)
		calc_label = __hash64String(inputfile)

		a = calc_label+'.in'
		new_inputfile = open(os.path.join(subdir,a),'w')
		new_inputfile.write(inputfile)
		new_inputfile.close()

		output_calcs[calc_label] = new_calcs[f]
		output_calcs[calc_label]['_AFLOWPI_INPUT_'] = inputfile
			
	except IOError as e:
		logging.error("%s not in %s" % (f,new_calcs[f]['_AFLOWPI_FOLDER_']))

    return output_calcs

#########################################################################################################################

def modifyCalcs(old,new,calcs):
	"""
	############
	##OBSOLETE##
	############
	Modify the calculation in each subdir and update the master dictionary
	
	Arguments:
	 old   (str) : string to replace
	 new   (str) : replacement string
	 calcs (dict): dictionary of dictionaries of calculations

	"""

	new_calcs = copy.deepcopy(calcs)
	output_calcs = collections.OrderedDict()
	
	for f,d in new_calcs.iteritems():
		try:
			subdir = new_calcs[f]['_AFLOWPI_FOLDER_']
			a = f+'.in'
			inputfile = open(os.path.join(subdir,a),'r').read()
			inputfile = re.sub(old,new,inputfile)
			calc_label = __hash64String(inputfile)
			
			a = calc_label+'.in'
			new_inputfile = open(os.path.join(subdir,a),'w')
			new_inputfile.write(inputfile)
			new_inputfile.close()

			output_calcs[calc_label] = new_calcs[f]
			output_calcs[calc_label]['_AFLOWPI_INPUT_'] = inputfile

		except IOError as e:
			logging.error("%s not in %s" % (f,new_calcs[f]['_AFLOWPI_FOLDER_']))

	return output_calcs

#####################################################################################################################



import atexit

def line_prepender(filename,new_text):
    ''' 
    prepends a file with a new line containing the contents of the string new_text.

    Arguments:
          filename (str): string of the filename that is to be prepended
	  new_text (str): a string that is one line long to be prepended to the file

    Returns:
          None

    '''

    try:
        with open(filename,'r+') as f:
            content = f.read()
            f.seek(0,0)
            f.write(new_text.rstrip('\r\n') + '\n' + content)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
############################################################################################################



############################################################################################################


def __findInBlock(oneCalc,ID,block=None,string=''):
    '''
    Looks for a string via a regular expression inside one of the command blocks in the _ID.py.
    Used to check before writing something as to avoid unintentionally having it written it twice.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          block (str): a string of the name of the command block in the _ID.py to look in.
	  string (str): the string to search for which can be in plain text or as a regex.

    Returns:
          True if the regex finds the pattern or False if it does not
          
    '''

    subdir = oneCalc['_AFLOWPI_FOLDER_']    
    fileName = '_'+ID+'.py'
    try:
        with open(os.path.join(subdir,fileName),'r') as inputfile:
            inputFileText = inputfile.read()

        try:
            isolatedString=re.split('#%s_BLOCK' % block,inputFileText)[-1]
            isolatedString=re.split('#END_%s_BLOCK' % block,inputFileText)[0]
        except:
            isolatedString=inptFileText

        return len(re.findall(string,isolatedString))!=0    
    
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        return False


def __calcsFromCalcList(calcList):
    '''
    ############
    ##OBSOLETE##
    ############
    Takes a list of dictionaries representing a calculation and turns the
    list of dictionaries into a dictionary of dictionaries

    Arguments:
          calcList (list): a list of dictionaries

    Returns:
          A dictionary of dictionaries
   
    '''
 
    largeSet=collections.OrderedDict()
    for calcs in calcList:
        calcCopy=copy.deepcopy(calcs)
        largeSet = dict(largeSet.items() + calcs.items())
    return largeSet


####################################################################################################################
###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ####
####################################################################################################################



# def useFunc(calcs,func=None,args=(),kwargs={},newStep=False):
    
#     sm=[]
#     badList=['__builtin__','__builtins__','_sh']
#     for var in dir(__main__):
#         try:
#             if inspect.ismodule(eval('__main__.%s'%var)):
#                 if var not in badList:
#                     sm.append('try:\n\timport %s\nexcept Exception,e:\n\tAFLOWpi.run._fancy_error_log(e)'%var)
#         except Exception,e:
#             AFLOWpi.run._fancy_error_log(e)

#     loadModString= '\n'.join(['%s'%x for x in sm])
    
#     funcLines=['\t%s'%x for x in inspect.getsourcelines(func)[0]]
# #    print funcLines[0].split(',')
#     print inspect.getargspec(func)
#     funcStr='try:\n'
#     funcStr+=''.join(funcLines)
#     funcStr+='\nexcept Exception,e:\n\tAFLOWpi.run._fancy_error_log(e)'

# #    print funcStr
#     callStr='try:\n\t'

#     callStr+='oneCalc = '
#     callStr+='%s('%func.__name__
#     argStr=','.join(['%s'%repr(x) for x in args])+','
#     kwargStr=','.join(['%s=%s'%(x[0],repr(x[1])) for x in kwargs.items()])
#     callStr+=argStr
#     callStr+=kwargStr
#     callStr+=')\nexcept Exception,e:\n\tAFLOWpi.run._fancy_error_log(e)'
# #    print callStr
#     print args

    
#     if newStep:
#         new_calcs = writeToScript(func.__name__,calcs)

#         __addToAll(new_calcs,block='IMPORT',addition=loadModString)
#         __addToAll(new_calcs,block='IMPORT',addition=funcStr)

#         return new_calcs

####################################################################################################################
###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ####
####################################################################################################################





#################################################################################################################














# def __loadAllCalcs(fileList):
#     newCalcs = collections.OrderedDict()
#     for item in fileList:
#         folder = os.path.dirname(item)
#         splitList=item.split('/')
#         calcID=splitList[-1][1:]

#         calcID=re.sub('.py','',calcID)

#         oneCalc = AFLOWpi.prep.__loadOneCalc(folder,calcID)
#         newCalcs[calcID]=oneCalc
    

#     return newCalcs


	##########################################################################################
	###NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED ###
	##########################################################################################








#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
######################################################################################################################### 




####################################################################################################################
###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ####
####################################################################################################################


############################################################################################################
#import AFLOWpi.pseudo as PT


####################################################################################################################
###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ###NOT FINISHED ####
####################################################################################################################


#     '''keep trying to save the calc to the calclog. if an error occurs because the file is already open by another job then wait 1 second and try again'''
#     while True:
#         try:
            
#             output = shelve.open(filename)
#             calcSave = copy.deepcopy(oneCalc)
#             output[ID]=calcSave
#             break 
#         except Exception,e:
#             AFLOWpi.run._fancy_error_log(e)
#             time.sleep(1)
#             pass

#     logging.debug('exiting __updatelos')




#

# def __swapPositions(symMatrix,order=[1,2,3]):
#     '''reorganize the atomic_positions to fit with the change in the cell vectors'''


#     temp=[]
#     tempPos=[]

    
#     '''transpose so we can switch the vertical vectors easier'''
#     transposedPositions=symMatrix.T
#     returnMatrix = copy.deepcopy(transposedPositions)

#     returnMatrix[0] = numpy.copy(transposedPositions[order[0]-1])
#     returnMatrix[1] = numpy.copy(transposedPositions[order[1]-1])
#     returnMatrix[2] = numpy.copy(transposedPositions[order[2]-1])
        
#     '''transpose the positions again so they're back to normal'''
#     return returnMatrix.T

# def __aflow2pw(inputString):


#     inputDict = AFLOWpi.retr.__splitInput(inputString)
#     try:
# 	    cell=AFLOWpi.retr.__cellStringToMatrix(inputDict['CELL_PARAMETERS']['__content__'])
	    
#     except:
# 	    '''if there's no cell_positions card don't worry about the convention'''
# 	    return inputString

#     positionMatrix = AFLOWpi.retr.__getPositions(inputString)
#     labels = AFLOWpi.retr.__getPosLabels(inputString)
#     pos,flag = AFLOWpi.retr.detachPosFlags(AFLOWpi.qe.regex.atomic_positions(inputString))


#     ibrav=AFLOWpi.retr.getIbravFromVectors(cell)

#     if ibrav==13:
#         CONV_MCLC=numpy.array([
#                 [1.0,  0.0,  0.0,],
#                 [0.0,  0.0,  1.0,],
#                 [0.0,  1.0,  0.0,],
#                 ])
#         cell=CONV_MCLC.dot(cell)
#         cell=AFLOWpi.retr.__prim2ConvVec(cell)



#         a = numpy.sqrt(cell[0].dot(cell[0].T)).getA()[0][0]
#         b = numpy.sqrt(cell[1].dot(cell[1].T)).getA()[0][0]
#         c = numpy.sqrt(cell[2].dot(cell[2].T)).getA()[0][0]
#         alpha = numpy.arccos(cell[1].dot(cell[2].T).getA()[0][0]/(b*c))*180.0/numpy.pi
#         beta  = numpy.arccos(cell[0].dot(cell[2].T).getA()[0][0]/(a*c))*180.0/numpy.pi
#         gamma = numpy.arccos(cell[0].dot(cell[1].T).getA()[0][0]/(a*b))*180.0/numpy.pi

# 	b_temp=b
# 	c_temp=c
# 	b=b_temp
# 	c=c_temp
	

# 	alpha_temp=alpha
# 	beta_temp=beta
# 	gamma_temp=gamma



# #	a=b_temp
# #	b=a_temp
# #	c=c_temp
# 	alpha=alpha_temp
# 	beta=beta_temp
# 	gamma=gamma_temp


	
#         MCLC_fixed=AFLOWpi.retr.abc2free(a,b,c,alpha=alpha,beta=beta,gamma=gamma,ibrav=ibrav)
#         MCLC_fixed=AFLOWpi.retr.__cellStringToMatrix(MCLC_fixed)

# 	positionMatrix=__swapPositions(positionMatrix,order=[1,3,2])


# 	pos=AFLOWpi.retr.__joinMatrixLabels(labels,positionMatrix)
# 	pos=AFLOWpi.retr.attachPosFlags(pos,flag)
# 	MCLC_fixed= AFLOWpi.retr.__cellMatrixToString(MCLC_fixed)
	



# 	inputDict['CELL_PARAMETERS']['__content__']=MCLC_fixed
# 	inputDict['ATOMIC_POSITIONS']['__content__']=pos

# 	inputString=AFLOWpi.retr.__joinInput(inputDict)

#     if int(ibrav)==4:
# 	    trans=numpy.array([
# 			    [ 1.,  1.,  0.,],
# 			    [-1.,  0.,  0.,],
# 			    [ 0.,  0.,  1.,],
# 			    ])
# 	    cell=trans.dot(cell)
# 	    cell= AFLOWpi.retr.__cellMatrixToString(cell)
# 	    inputDict['CELL_PARAMETERS']['__content__']=cell
# 	    inputString=AFLOWpi.retr.__joinInput(inputDict)

    

#     return inputString


