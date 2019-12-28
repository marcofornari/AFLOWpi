import os
import datetime
import pickle
import logging 
import re 
import subprocess
import sys                                                                     
import signal
import collections
import copy
import logging.handlers
import AFLOWpi
import time
import queue
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
import tempfile
import pipes


if __name__!='__main__':
    engineDict={}

def _exitClean(signal,frame):
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

def _recordDeath(signal,frame):
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

#    index=globals()["__TEMP__INDEX__COUNTER__"]
    try:
#        if int(__main__.oneCalc['_AFLOWPI_INDEX_'])==index:
            __main__.oneCalc['__status__']['Error']="Killed"
            AFLOWpi.prep._saveOneCalc(__main__.oneCalc,__main__.ID)
    except:
        pass
    sys.exit(0)
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

    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.prep._addToBlock(oneCalc,ID,'CLEANUP','\nimport os\nos.system("rm -rf %s/_*")\n' % oneCalc['_AFLOWPI_FOLDER_'])            
        except Exception as e:
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
    
    if '__atexitList' not in list(globals().keys()):
        #######################################
        @atexit.register
        def _runatexitList():
            for item in __atexitList:
                item[0](*item[1],**item[2])
        global __atexitList
        __atexitList = []

    __atexitList.append((command,args,kwargs))




def _colorize_message(string,level='ERROR',show_level=True):
    """
    Colorizes text. Used for colorizing the logging
    
    Arguments:
          string (str): A string of text

    Keyword Arguments:
          level (str): Specific colors are chosen for logging message type

    Returns:
          levelname_color (str): Colorized version of the string input
          
    """

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = list(range(8))

    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[1;%dm"
    BOLD_SEQ = "\033[1m"

    COLORS = {
        'WARNING': YELLOW,
        'INFO': WHITE,
        'DEBUG': BLUE,
        'CRITICAL': YELLOW,
        'ERROR': RED,
        'GREEN': GREEN        ,
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

    logging.error(AFLOWpi.run._colorize_message('%%%DEBUG%%%'))
    logging.error(e)
    _, _, tb = sys.exc_info()
    errorList =  traceback.format_list(traceback.extract_tb(tb)[-6:])[-6:]
    for errorMSG in errorList:
        logging.error(AFLOWpi.run._colorize_message(errorMSG))     
        logging.error(AFLOWpi.run._colorize_message('%%%DEBUG%%%'))       

    try:
        if AFLOWpi.prep._ConfigSectionMap('prep','log_level').upper() == 'DEBUG':
            print('%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%')
            try:
                print(e)
                for errorMSG in errorList:
                    print(errorMSG)

            except Exception as e:
                print(e)

            print('%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%')
    except:
        pass


def _getEnginePath(engine,calcType):
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
        execPath=''
    return execPath

################################################################################################################

################################################################################################################
def _getExecutable(engine,calcType):
    """
    Gives the name of the executable file for the ab initio engine
    for a given type of calculation.

    OBSOLETE. NEEDS REMOVAL. AFLOWpi.run._getEnginePath is almost identical 
    
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
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        execPath=''
    return executable



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

    calcName=AFLOWpi.run._get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])

    stderrFileNm=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'_cluster'+'.stderr')
    stdoutFileNm=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'_cluster'+'.stdout')

    qsubRef = re.sub('_AFLOWPI_CALCNM_',calcName,qsubRefFileString)
    qsubRef = re.sub('_AFLOWPI_STDOUT_',stdoutFileNm,qsubRef)
    qsubRef = re.sub('_AFLOWPI_STDERR_',stderrFileNm,qsubRef)

    clusterTypeDict['qsubRef']=qsubRef
    logging.debug('exiting generateSubRef')


    return clusterTypeDict



#############################################################################################################

#############################################################################################################
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

    print('EMR NOT IMPLEMENTED')
    raise SystemExit


    testOne(calcs,calcType='emr',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    gipawdir=AFLOWpi.prep._ConfigSectionMap('prep','gipaw_dir')
    if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print(('ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir))
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)

    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='emr')
        except Exception as e:
            _fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='emr')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


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
    print('HYPERFINE NOT IMPLEMENTED')
    raise SystemExit


    testOne(calcs,calcType='hyperfine',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)

    gipawdir=AFLOWpi.prep._ConfigSectionMap('prep','gipaw_dir')
    if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print(('ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir))
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)


    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='hyperfine')
        except Exception as e:
            _fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='hyperfine')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)



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

    print('NMR NOT IMPLEMENTED')
    raise SystemExit

    testOne(calcs,calcType='nmr',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)
    gipawdir=AFLOWpi.prep._ConfigSectionMap('prep','gipaw_dir')
    if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print(('ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir))
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)



    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='nmr')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='nmr')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)




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

    print('GVECTORS NOT IMPLEMENTED')
    raise SystemExit

    testOne(calcs,calcType='gvectors',engine=engine,execPrefix=execPrefix,execPostfix=execPostfix,holdFlag=holdFlag,config=config)

    gipawdir=AFLOWpi.prep._ConfigSectionMap('prep','gipaw_dir')
    if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
        symlink=True
    else:
        symlink=False
    gipawExLoc = os.path.join(gipawdir,'gipaw.x')
    if not os.path.exists(gipawExLoc):
        logging.error('ERROR: gipaw not found in %s please check your config file. EXITING' % gipawdir)
        print(('ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir))
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)


    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='gvectors')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='gvectors')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


#############################################################################################################
#############################################################################################################
def _skeletonRun(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None):
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

def scf(calcs,engine='',execPrefix=None,execPostfix=None,holdFlag=True,config=None,exit_on_error=True):
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
    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',exit_on_error=exit_on_error)
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


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
    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='dos')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='dos')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)

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
    for ID,oneCalc in list(calcs.items()):
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='pdos')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='pdos')
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)



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
        engine = AFLOWpi.prep._ConfigSectionMap("run",'engine')
    
    testOne(calcs,calcType='bands',engine=engine,execPrefix=execPrefix,execPostfix='',holdFlag=holdFlag,config=config)
    for ID,oneCalc in list(calcs.items()):
#        try:
#            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='bands')
#        except Exception as e:
#            AFLOWpi.run._fancy_error_log(e)
#        try:
#            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='bands')
#        except Exception as e:
#            AFLOWpi.run._fancy_error_log(e)

        bands_pp_run_string = '''if oneCalc['__execCounter__'] <=%s:
        AFLOWpi.run._bands_pp(__submitNodeName__,oneCalc,ID)
        oneCalc['__execCounter__']+=1
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
'''%oneCalc['__execCounterBkgrd__']
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',bands_pp_run_string)
        oneCalc['__execCounterBkgrd__']+=1

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
                if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
                    execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

                else:
                    execPrefix=''

            if execPostfix == None:
                if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
                    execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")

                else:
                    execPostfix=''
            
            if calcType=='bands':
                if len(re.findall(r'pool',execPostfix))!=0:
                    execPostfix=' '

            if engine=='':
                engine = AFLOWpi.prep._ConfigSectionMap("run",'engine')

            engineDir  = AFLOWpi.prep._ConfigSectionMap("prep",'engine_dir')	
            if os.path.isabs(engineDir) == False:
                configFileLocation = AFLOWpi.prep._getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                engineDir =  os.path.join(configFileLocation, engineDir)



            try:
                enginePath =  AFLOWpi.run._getEnginePath(engine,calcType)
            except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)

            for files in enginePath:

                if os.path.exists(os.path.join(engineDir,files))==False:
                    print(('%s not found. Check your config file to make sure engine_dir path is correct and that %s is in that directory..Exiting' %(files,files)))
                    logging.error('%s not found. Check your config file to make sure engine_dir path is correct and that %s is in that directory..Exiting' %(files,files))
                    raise SystemExit

                if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
                    AFLOWpi.prep.totree(os.path.abspath(os.path.join(engineDir,files)),calcs)
                else:
                    for ID,oneCalc in list(calcs.items()):
                        try:
                            os.symlink(os.path.abspath(os.path.join(engineDir,files)),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],files))
                        except OSError:
                            os.unlink(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],files))
                            os.symlink(os.path.abspath(os.path.join(engineDir,files)),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],files))

            logfile = os.path.abspath(os.path.join((os.path.dirname(calcs[random.choice(list(calcs.keys()))]['_AFLOWPI_FOLDER_'])),'AFLOWpi','LOG.log'))
        except Exception as e:
            _fancy_error_log(e)
        
        #############################################################################################################
        try:
            calcsCopy=copy.deepcopy(calcs)
            newCalcs = collections.OrderedDict()

            for ID,oneCalc in list(calcs.items()):
                CONFIG_FILE='../AFLOWpi/CONFIG.config'
                calcs[ID]['_AFLOWPI_CONFIG_']=CONFIG_FILE
                calcs[ID]['calcType']=calcType                

        except Exception as e:
            _fancy_error_log(e)
        #############################################################################################################
        clusterType = AFLOWpi.prep._ConfigSectionMap("cluster",'type')
        if clusterType=='':
            clusterType=''
        try:
            if 'firstCalcList' not in list(globals().keys()) or holdFlag==False:

                try:
                    baseID = ID.split('_')[0]
                    try:
                        splitIndexInt = int(oneCalc['__TEMP__INDEX__COUNTER__'])
                        prev_ID = '%s_%.02d' % (baseID.split('_')[0],splitIndexInt-1)
                    except:
                        prev_ID=ID
                        pass                   

                    if os.path.exists(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+prev_ID+'.py')):
                        pass

                    else:
                        prev_ID=ID

                except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)

            elif firstCalcList!=False:
                pass


            firstCalcList=False
        except Exception as e:
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

    for ID,oneCalc in list(calcs.items()):
        if oneCalc['__chain_index__']!=1:
            log_file=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.oneCalc' % ID)            
            os.remove(log_file)
        else:
            oneCalc['__execCounter__']=0
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)



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


def _get_qsub_name(path):
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

    try:
        calcName = '_'.join([j for j in c_name.split('_')[-4:] if j!=''])
    except Exception as e:
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


#

def _qsubGen(oneCalc,ID):
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

    logging.debug('ENTERING _qsubGen')
    calcName=AFLOWpi.run._get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])
    execFileString = os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py'))
    tmpdir_envar=''
    ls_option = AFLOWpi.prep._ConfigSectionMap("cluster",'local_scratch').strip().lower()

    if ls_option=='true':
        tmpdir = AFLOWpi.prep._ConfigSectionMap("cluster",'local_scratch_dir').strip()
        if tmpdir!='':
            tmpdir_envar="export TMPDIR='%s'\n"%tmpdir

    try:
        qsubFilename = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.qsub' % ID)
        prevQsubString = ''
        clusterType = AFLOWpi.prep._ConfigSectionMap("cluster",'type').strip().upper()

        if clusterType in ['PBS','UGE']:
            qsubRefFileName = AFLOWpi.prep._ConfigSectionMap("cluster",'job_template')
            if os.path.isabs(qsubRefFileName) == False:
                configFileLocation = AFLOWpi.prep._getConfigFile()
                configFileLocation = os.path.dirname(configFileLocation)
                qsubRefFileName =  os.path.join(configFileLocation, qsubRefFileName)        
            qsubRefString = file(qsubRefFileName,'r').read()
            qsubRef = generateSubRef(qsubRefString,oneCalc,ID)['qsubRef']
            qsubRef = AFLOWpi.prep.remove_blank_lines(qsubRef)



            #put the name of the AFLOWpi job in the  PBS submit script
            if clusterType=='UGE':
                name_regex=re.compile(r'\s*#$.*-[nN].*\n')
                calcNameString = '\n#$ -N %s\n'%calcName
            else:
                name_regex=re.compile(r'\s*#PBS.*-[nN].*\n')
                calcNameString = '\n#PBS -N %s\n'%calcName

            if len(name_regex.findall(qsubRef)):
                qsubRef=name_regex.sub(calcNameString,qsubRef)
            else:
                qsubRef=calcNameString+qsubRef


        #not sure if this will work on UGE so we'll leave it only for PBS
        if clusterType=='PBS':
                dash_e_regex=re.compile(r'#PBS\s+-e\s+.*\n')
                dash_o_regex=re.compile(r'#PBS\s+-o\s+.*\n')
                dash_oe_regex=re.compile(r'#PBS\s+-oe\s+.*\n')
                stderr_name=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'_cluster'+'.stderr')
                stdout_name=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'_cluster'+'.stdout')
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




        with open(qsubFilename,'w') as qsubFile:
            qsubSub='''%scd %s
python %s''' % (tmpdir_envar,oneCalc['_AFLOWPI_FOLDER_'],os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py'))
            #redirect the -e, -o or -oe to the directory for that job
            if len(re.findall('_AFLOWPI_QSUB_',qsubRef)):
                qsubRef = re.sub('_AFLOWPI_QSUB_',qsubSub,qsubRef)
            else:
                qsubRef+='\n'+qsubSub+'\n'

            qsubFileString=qsubRef
            qsubFile.write(qsubFileString)
            return qsubFileString



    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        raise SystemExit



#############################################################################################################
def _testOne(ID,oneCalc,engine='',calcType='',execPrefix=None,execPostfix=None): 
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
        if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
            execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")


        else:
            execPrefix=''

    if execPostfix == None:
        if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
            execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
        else:
            execPostfix=''


    logging.debug('entering _testOne') 
    try:
        global configFile
        configFile = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'../AFLOWpi/CONFIG.config')

        AFLOWpi.prep._forceGlobalConfigFile(configFile)

        if execPrefix == None:
            if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
                execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
            else:
                execPrefix=''

        if execPostfix == None:
            if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
                execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
            else:
                execPostfix=''
    except Exception as e:
        _fancy_error_log(e)

    try:
        engineDir  = AFLOWpi.prep._ConfigSectionMap("prep",'engine_dir')
        if os.path.isabs(engineDir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            engineDir =  os.path.join(configFileLocation, engineDir)
    except Exception as e:
        pass


    ri = ID+'.in'
    ro = ID+'.out' 

    rd = oneCalc['_AFLOWPI_FOLDER_']
    configFile = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'../AFLOWpi/CONFIG.config')

    logging.debug('exiting _testOne')

#############################################################################################################


def _submitJob(ID,oneCalc,__submitNodeName__,sajOverride=False,forceOneJob=False):
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

    logging.debug('entering _submitJob')
    folder = oneCalc['_AFLOWPI_FOLDER_']

    try:
        clusterType = AFLOWpi.prep._ConfigSectionMap("cluster","type").upper()

        if len(clusterType)==0:
            clusterType='None'
        elif clusterType==None:
            clusterType='None'
    except Exception as e:
        _fancy_error_log(e)

    try:
        stepsAsJobs = AFLOWpi.prep._ConfigSectionMap("cluster","steps_as_jobs").lower()
    except Exception as e:
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
        print('Submitting steps as one cluster job per calculation chain.')
        logging.info('Submitting steps as one cluster job per calculation chain.')

    if (socket.gethostname() != __submitNodeName__ and stepsAsJobs.lower()!='true') or forceOneJob==True:
        if sajOverride == False:
            clusterType=''

###########################################################################################

    if ID == oneCalc['_AFLOWPI_PREFIX_'][1:]:
        pass
    else:
        ID='_'.join(ID.split('_')[:3])

    copyback_level = AFLOWpi.prep._ConfigSectionMap("cluster","copyback_every_step").lower()
    if copyback_level.strip().lower()=='true':
        AFLOWpi.prep._from_local_scratch(oneCalc,ID)

    try:
        #for calc subset...
        if oneCalc['_AFLOWPI_FOLDER_']==__main__.oneCalc['_AFLOWPI_FOLDER_']:
            #submit the cluster submission 
            if os.path.exists(os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.qsub'))):
                submit_ID=ID
            else:
                submit_ID=__main__.ID
        else:
            submit_ID=ID
    except Exception as e:
        submit_ID=ID

    submit_ID=ID

    cluster_daemon=False
    cluster_daemon_option = AFLOWpi.prep._ConfigSectionMap("cluster","daemon").lower()
    if cluster_daemon_option=='true':
        cluster_daemon=True
        
    if clusterType.upper() in ['PBS','SLURM','UGE']:
        '''if we're using a daemon then write the .submit file and return'''
        if cluster_daemon==True:
            cluster_submit_file = os.path.abspath(os.path.join(folder,'_'+submit_ID+'.submit'))            
            with open(cluster_submit_file,'w') as cluster_submit_file_obj:
                cluster_submit_file_obj.write('%s'%submit_ID)
            return

        #if we're not copying back every step we need to copy back here for a
        #new cluster job submission to make sure we have the files to restart
        if copyback_level.strip().lower()!='true':
            AFLOWpi.prep._from_local_scratch(oneCalc,ID)

        additional_flags="" 
#        if clusterType.upper() in ['PBS','SLURM']:
#            additional_flags="" 
        
            
        try:
            try:
                submitCommand = AFLOWpi.run._get_cluster_submit_command()
                
            except Exception as e:
                _fancy_error_log(e)
            
            qsubFilename = os.path.abspath(os.path.join(folder,'_'+submit_ID+'.qsub'))

            calcName=AFLOWpi.run._get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])

            queue = ''
            if AFLOWpi.prep._ConfigSectionMap("cluster","queue") !='':
                queue = '-q %s ' %  AFLOWpi.prep._ConfigSectionMap("cluster","queue")
        except Exception as e:
            _fancy_error_log(e)
        try:
            
            command  = "ssh -o StrictHostKeyChecking=no %s '%s %s -N %s %s %s' " % (__submitNodeName__,submitCommand,queue,calcName,qsubFilename,additional_flags)
            job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True)
            job.communicate()               

            logging.info('submitted %s to the queue' % submit_ID)
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
            logging.info("Trying to submit job with no -N option")

            try:
                command = "ssh -o StrictHostKeyChecking=no %s '%s %s %s' " % (__submitNodeName__,submitCommand,qsubFilename,additional_flags)
                job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True)
                job.communicate()               

            except Exception as e:
                logging.info("Trying to submit job without ssh and -N option")
                try:
                    baseID = ID.split('_')[0]
                    baseID = ID
                    qsubFilename = os.path.abspath(os.path.join(folder,'_'+baseID+'.qsub'))
                    os.system("ssh -o StrictHostKeyChecking=no %s '%s %s -N %s %s %s' " % (__submitNodeName__,submitCommand,queue,calcName,qsubFilename,additional_flags))
                except Exception as e:
                    logging.info("Trying to submit job without ssh and -N option")
                    try:
                        baseID = ID.split('_')[0]
                        baseID = ID
                        qsubFilename = os.path.abspath(os.path.join(folder,'_'+baseID+'.qsub'))
                        os.system("%s %s %s %s" % (submitCommand,queue,qsubFilename,additional_flags))
                    except Exception as e:
                        logging.info("Trying to submit job without ssh and -N option")
                        try:
                            os.system("%s %s %s %s" % (submitCommand,queue,qsubFilename,additional_flags))
                        except Exception as e:
                            _fancy_error_log(e)

    else:
        try:

            execFileString  = os.path.abspath(os.path.join(folder,'_'+ID+'.py'))
            
            exec(compile(open(execFileString, "rb").read(), execFileString, 'exec'))

        except KeyboardInterrupt:
            print('Got Keyboard Exit signal, exiting')
            logging.debug('Got Keyboard Exit signal, exiting')
            sys.exit(0)
        except Exception as e:
            _fancy_error_log(e)

    logging.debug('exiting _submitJob')


##############################################################################################

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

    cluster_daemon_option = AFLOWpi.prep._ConfigSectionMap("cluster","daemon").lower()
    if cluster_daemon_option=='true':
        AFLOWpi.run._generate_submission_daemon_script(calcs)

    logging.debug('entering submitFirstCalcs__')
    if '__submitNodeName__' not in list(globals().keys()):
        global __submitNodeName__
        __submitNodeName__ = socket.gethostname()
    try:
        __submitNodeName__=__main__.__submitNodeName__
    except:
        pass

    logging.debug('sending from %s' % __submitNodeName__)
    
    if AFLOWpi.prep._ConfigSectionMap('cluster','type') != '':
        submit_message = '\n*** Submitting Workflow ***\n'
    else:
        submit_message = '\n*** Starting Workflow ***\n'

    print((AFLOWpi.run._colorize_message(submit_message,level='ERROR',show_level=False)))

    for ID,oneCalc in list(calcs.items()):
        AFLOWpi.run._submitJob(ID,oneCalc,__submitNodeName__)
    logging.debug('exiting submitFirstCalcs__')        


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
def _restartPW(oneCalc,ID,a,__submitNodeName__,):
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
    logging.debug('Entering AFLOWpi.run._restartPW')

    try:
        walltime = AFLOWpi.run._grabWalltime(oneCalc,ID)
    except:
        return 

    with open('%s.in' % ID,'r') as inputFileObj:
        inFileString= inputFileObj.read()
    inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    #bring wfc,hub, or any other files needed from local scratch if need be
    walltimeSec=90000000000
    try:
        walltimeSec = int(inputDict["&control"]["max_seconds"])
    except:
        '''if we can't find max_seconds in the input we just exit.'''
        logging.debug('Not pwscf. Exiting AFLOWpi.run._restartPW')
        return

    try:
        """see if calc reached it's end."""
        if (time.time()-a)>walltimeSec:
            inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
            """if it did then set to restart the calc"""

            inputDict["&control"]["restart_mode"] = "'restart'"
            ''' if the calc never had the time to start then just restart it from scratch'''

            inputDict["&control"]["max_seconds"] = int(walltimeSec)

            restartInput = AFLOWpi.retr._joinInput(inputDict)
            '''write the new input file with the restarting'''
            with open('%s.in' % ID,"w") as inputfile:
                inputfile.write(restartInput)                  
            '''save to _<ID>.oneCalc'''
            oneCalc['_AFLOWPI_INPUT_']=restartInput
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            """now resubmit"""
            resubmitID = oneCalc['__qsubFileName__'][1:-5]

            #resubmit
            '''if we're at within 15 seconds of the end of walltime. don't submit and let restartScript do it'''

            oneCalc['__status__']['Restart']+=1

            AFLOWpi.run._submitJob(resubmitID,oneCalc,__submitNodeName__,sajOverride=True)                             
            logging.debug('Restart Job Submitted. Exiting AFLOWpi.run._restartPW')
            time.sleep(10)
            '''and exit cleanly.'''
            sys.exit(0)

        else:
            logging.debug('Job completed in time. No need for restart. Exiting AFLOWpi.run._restartPW')

    except Exception as e:
        _fancy_error_log(e)


def _restartGIPAW(oneCalc,ID,a,__submitNodeName__,):
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
    inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    try:
        walltimeSec = int(inputDict["&control"]["max_seconds"])
    except:
        '''if we can't find max_seconds in the input we just exit.'''
        sys.exit(0)
    try:
        """see if calc reached it's end."""
        if (time.time()-a)>walltimeSec:
            inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
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
            restartInput = AFLOWpi.retr._joinInput(inputDict)
            '''write the new input file with the restarting'''
            with open('%s.in' % ID,"w") as inputfile:
                inputfile.write(restartInput)                  
            '''save to _<ID>.oneCalc'''
            oneCalc['_AFLOWPI_INPUT_']=restartInput
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            """now resubmit"""
            resubmitID = oneCalc['__qsubFileName__'][1:-5]

            oneCalc['__status__']['Restart']+=1
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)

            AFLOWpi.run._submitJob(resubmitID,oneCalc,__submitNodeName__,sajOverride=True)                 
            '''and exit cleanly.'''
            sys.exit(0)
    except Exception as e:
        _fancy_error_log(e)



################################################################################################################
def _getWalltime(oneCalc,ID):
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

    cluster_type=AFLOWpi.prep._ConfigSectionMap("cluster","type").upper()

    try:
        with open(oneCalc['__qsubFileName__'],"r") as qsubfileObj:
            qsubFileString = qsubfileObj.read()
    except Exception as e:
        return 50000000
    try:
        if cluster_type=='PBS':
            walltime = re.findall("walltime\s*=\s*([0-9:]*)",qsubFileString)[-1]

        elif cluster_type=='SLURM':
            walltime = re.findall("--time\s*=\s*([0-9:]*)",qsubFileString)[-1]
        else:
            walltime='3000:00:00'
    except:
        walltime='3000:00:00'

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


def _generic_restart_check(oneCalc,ID,__submitNodeName__):
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
        walltime,startTime=AFLOWpi.retr._grabWalltime(oneCalc,ID)    
    except:
        return

    percent_timer=0.90
    from_config = AFLOWpi.prep._ConfigSectionMap("cluster","restart_buffer").lower()
    if from_config!='':
        percent_timer = float(from_config)
    
    if percent_timer<1.0:
        num_secs=int(float(walltime)*percent_timer)
    else:
        num_secs=int(float(walltime)-percent_timer)

    try:
        timediff=time.time()-startTime
    except Exception as e:
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
            
    resubmitID = oneCalc['__qsubFileName__'][1:-5]

    if (runTime-60.0)<0.0:
        ID = __get_index_from_pp_step(oneCalc,ID)
        '''if there's not enough time to do anything so resubmit before anything runs'''
        '''make sure the scratch is copied back'''
        logging.info('job only has 60 seconds left. attempting to copy from local scratch if needed and resubmit job')
        '''force a new job'''
        AFLOWpi.run._submitJob(resubmitID,oneCalc,__submitNodeName__,sajOverride=True)                 
        logging.debug('Not enough time to start next step. Restarting')
        time.sleep(5)
        sys.exit(0)        



def _setupRestartPW(oneCalc,ID,__submitNodeName__):
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
        walltime,startTime=AFLOWpi.run._grabWalltime(oneCalc,ID)
    except:
        '''if no walltime log in oneCalc (running local mode)'''
        return oneCalc,ID

    '''set buffer to time before walltime to stop the job'''
    percent_timer=0.90
    from_config = AFLOWpi.prep._ConfigSectionMap("cluster","restart_buffer").lower()
    if from_config!='':
        percent_timer = float(from_config)


        
    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"r") as inputfile:
        inputStr=inputfile.read()
    '''try to split input to include max time if it doesn't work stop and return'''
    try:
        inputDict = AFLOWpi.retr._splitInput(inputStr)
        calc_type=inputDict["&control"]["calculation"].strip(""" '`" """)
    except:
        return oneCalc,ID
    '''ibef calc type was found (i.e. not post processing) try to include max_seconds'''
    typeList = ['bands','nscf','scf','vc-relax','relax']

    try:
        timediff=time.time()-startTime
    except Exception as e:
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
                    resubmitID = oneCalc['__qsubFileName__'][1:-5]
                    ID = __get_index_from_pp_step(oneCalc,ID)
                    '''force a new job'''
                    AFLOWpi.run._submitJob(resubmitID,oneCalc,__submitNodeName__,sajOverride=True)                 
                    logging.debug('Not enough time to start pwscf. Restarting')
                    time.sleep(10)
                    sys.exit(0)        
            
            inputDict["&control"]["max_seconds"] = "%s" % int(runTime)
            restartInput = AFLOWpi.retr._joinInput(inputDict)
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"w") as inputfile:
                inputfile.write(restartInput)
            oneCalc['_AFLOWPI_INPUT_']=restartInput
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)

        except Exception as e:
            logging.debug('if option to restart doesnt exist in the input of %s this will be thrown' % ID)
            _fancy_error_log(e)

        return oneCalc,ID



def _setupRestartGIPAW(oneCalc,ID):
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
        walltime,startTime=AFLOWpi.run._grabWalltime(oneCalc,ID)
    except:
        return oneCalc,ID
    
    globals()['__GLOBAL_CLOCK__']=startTime

    percent_timer=0.90
    try:
        percent_timer = AFLOWpi.prep._ConfigSectionMap("cluster","restart_buffer").lower()
        percent_timer = float(percent_timer.strip())
    except:
        percent_timer=0.90

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"r") as inputfile:
        inputStr=inputfile.read()
    '''try to split input to include max time if it doesn't work stop and return'''
    try:
        inputDict = AFLOWpi.retr._splitInput(inputStr)
        calc_type=inputDict["&inputgipaw"]["calculation"].strip(""" '`" """)
    except:
        return oneCalc,ID
    '''ibef calc type was found (i.e. not post processing) try to include max_seconds'''
    typeList = ['emr','efg','g_vectors','hyperfine']
    if calc_type in typeList:
        try:
            timediff=time.time()-globals()['__GLOBAL_CLOCK__']

        except Exception as e:
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
            restartInput = AFLOWpi.retr._joinInput(inputDict)
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"w") as inputfile:
                inputfile.write(restartInput)
            oneCalc['_AFLOWPI_INPUT_']=restartInput
        except Exception as e:
            logging.debug('if option to restart doesnt exist in the input of %s this will be thrown' % ID)
            _fancy_error_log(e)

        return oneCalc,ID




def _writeWalltimeLog(oneCalc,ID,walltimeDict):
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
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
    except:
        pass

def _readWalltimeLog(oneCalc,ID):
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



def _grabWalltime(oneCalc,ID):
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
        output = AFLOWpi.run._readWalltimeLog(oneCalc,ID)

        AFLOWpi.run._writeWalltimeLog(oneCalc,ID,output)

        '''grab the first walltime in case stepsasjobs==false for all jobs in chain'''
        walltimeStart=output['walltime']    
        startTime=output['start']

        logging.debug('startTime Test: %s' %startTime)
    except Exception as e:
        logging.debug(e)
        raise SystemExit
        return float(40000000),float(0)
    return float(walltimeStart),float(startTime)

def _setStartTime(oneCalc,ID,__CLOCK__):
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
    
#    if AFLOWpi.prep._ConfigSectionMap("cluster","type").lower()=='':
#        return oneCalc
    try:
        try:
            walltime_dict=oneCalc['__walltime_dict__']
        except:
            oneCalc['__walltime_dict__']={}
            walltime_dict={}

        if  '__GLOBAL_CLOCK__' not in list(globals().keys()):
            '''only move from local scratch if this script is starting fresh'''

            AFLOWpi.prep._to_local_scratch(oneCalc,ID)
            oneCalc['__walltime_dict__']={}
            globals()['__GLOBAL_CLOCK__']=__CLOCK__
            walltime_dict['walltime'] = AFLOWpi.run._getWalltime(oneCalc,ID)         
            walltime_dict['start']=__CLOCK__
            AFLOWpi.run._writeWalltimeLog(oneCalc,ID,walltime_dict)
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        try: 

            oneCalc['__walltime_dict__']={}
            walltime_dict={}
            walltime_dict['walltime'] = 40000000000.0         
            walltime_dict['start']=0.0

            AFLOWpi.run._writeWalltimeLog(oneCalc,ID,walltime_dict)
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)

    return oneCalc


def _restartScript(oneCalc,ID,PID):
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
        percent_timer = AFLOWpi.prep._ConfigSectionMap("cluster","restart_buffer").lower()
        percent_timer = float(percent_timer.strip())
    except:
        percent_timer=0.90

    walltime,startTime = AFLOWpi.run._grabWalltime(oneCalc,ID)
    globals()['__GLOBLAL_CLOCK__']=startTime
    
    walltime  =  int(walltime)
    
    __submitNodeName__=__main__.__submitNodeName__

    stepsasjobs = AFLOWpi.prep._ConfigSectionMap("cluster","steps_as_jobs").lower()
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
        except Exception as e:
            logging.debug(e)
        time.sleep(1.0)

        try:
            '''check to see if main script is running..if not it's 
           not then something went wrong in the main script and this should also exit'''

            os.kill(PID, 0) 
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
            logging.info('main script did not exit cleanly. killing __restartScript subprocess')
            print('main script did not exit cleanly. killing __restartScript subprocess')
            sys.exit(0)


    logging.info('Job not completed in time. Restarting from that step.')

    oneCalc =AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],ID)
    try:
        oneCalc['__status__']['Restart']+=1
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
        AFLOWpi.run._submitJob(ID,oneCalc,__submitNodeName__,sajOverride=True)
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
    logging.info('Resubmission successful. Waiting 5 seconds then ending main script')
    '''after the job is submitted wait ten seconds to exit the script'''
    time.sleep(60)
    '''and exit gracefully'''
    logging.info('telling main process to stop')
    os.kill(PID,signal.SIGUSR1) 
    sys.exit(0)


##################################################################################################################
def _convert_fortran_double(fort_double_string,string_output=False):
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
            print('could not convert fortran float string to numpy float')
            logging.info('could not convert fortran float string to numpy float')
            return fort_double_string
    if len(fort_double_list)>2:
        print('could not convert fortran float string to numpy float')
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

    except Exception as e:
        print(e)
        print('could not convert fortran float string to numpy float')
        logging.info('could not convert fortran float string to numpy float')
        return fort_double_string




def _get_index_from_pp_step(oneCalc,ID):
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
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        return ID


def _PW_bands_fix(oneCalc,ID):
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
    except Exception as e:
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
                except Exception as e:
                    pass
                total_path_length=numpy.asarray([float(x) for x in i.split()])
            else:
                try:
                    k_val.extend(numpy.asarray([float(x) for x in i.split()]))
                except Exception as e:
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

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)




def _vcrelax_error_restart(ID,oneCalc,__submitNodeName__):
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


    submit_ID = AFLOWpi.run._get_index_from_pp_step(oneCalc,ID)
    try:        
        if AFLOWpi.prep._check_restart(oneCalc,ID):
            oneCalc,submit_ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','restart_mode',"'from_scratch'")
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
                        oneCalcBase =AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],submit_ID)
                        oneCalcBase['__runList__']=[]
                        oneCalc['__runList__']=[]
                        AFLOWpi.prep._saveOneCalc(oneCalcBase,submit_ID)
                except:
                    pass

                AFLOWpi.run._submitJob(submit_ID,oneCalc,__submitNodeName__,forceOneJob=True)
                

def _io_error_restart(ID,oneCalc,__submitNodeName__):
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
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','restart_mode',"'from_scratch'")
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    error_regex_list=[]
    error_regex_list.append(re.compile(r'Error in routine diropn'))
    error_regex_list.append(re.compile(r'Error in routine davcio'))
    error_regex_list.append(re.compile(r'Error in routine pw_readfile'))
    error_regex_list.append(re.compile(r'Error in routine seqopn'))
    error_regex_list.append(re.compile(r'kpt grid not Monkhorst-Pack'))
    submit_ID = AFLOWpi.run._get_index_from_pp_step(oneCalc,ID)
    with open('%s.out'% ID,'r') as outFileObj:
        outFileString = outFileObj.read()
        for re_error in  error_regex_list:
            if len(re_error.findall(outFileString))!=0:
                if re_error.findall(outFileString)=='Error in routine pw_readfile' or re_error.findall(outFileString)=='kpt grid not Monkhorst-Pack':
                    try:
                        oneCalcBase =AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],submit_ID)
                        oneCalcBase['__runList__']=[]
                        oneCalc['__runList__']=[]
                        AFLOWpi.prep._saveOneCalc(oneCalcBase,submit_ID)
                    except:
                        pass
                try:

                    logging.warning("%s: Trying again..."%re_error.findall(outFileString)[-1])
                except:
                    pass
                
                try:
                    oneCalcBase =AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],submit_ID)

                    restart_num = oneCalcBase['__status__']['Restart']
                    oneCalcBase['__status__']['Restart']=restart_num+1
                    oneCalc['__status__']['Restart']=restart_num+1
                    if restart_num>10:
                        oneCalcBase["__execCounter__"]=0
                        oneCalc["__execCounter__"]=0
                    if restart_num>15:
                        sys.exit(0)

                    AFLOWpi.prep._saveOneCalc(oneCalcBase,submit_ID)
                except:
                    pass

                AFLOWpi.run._submitJob(submit_ID,oneCalc,__submitNodeName__,forceOneJob=True)
                sys.exit(0)

##################################################################################################################

def _qe__pre_run(oneCalc,ID,calcType,__submitNodeName__,engine):
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

    clusterType = AFLOWpi.prep._ConfigSectionMap("cluster","type")
    home = os.path.abspath(os.curdir)
    rd ='./'

    if calcType=='scf':    
        oneCalc,ID = AFLOWpi.prep._setup_local_scratch(oneCalc,ID)

    try:
        if clusterType.upper() in ['PBS','SLURM','UGE']:
            '''if there's 60 seconds or less left exit immediately and don't start a new job'''

            AFLOWpi.run._generic_restart_check(oneCalc,ID,__submitNodeName__)
            if calcType in ['scf']:    
                oneCalc,ID = AFLOWpi.run._setupRestartPW(oneCalc,ID,__submitNodeName__)



            if calcType in ['emr','nmr','hyperfine','gvectors']:
                oneCalc,ID,restartBool = __setupRestartGIPAW(oneCalc,ID)




    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    try:
        if calcType!='scf' and calcType!=None and calcType!='custom':
                with open(os.path.join(rd,"%s_%s.in") % (ID,calcType),'w+') as PPInput:
                    PPInput.write(AFLOWpi.run._makeInput(oneCalc,engine,calcType,ID=ID))

    except Exception as e:
        _fancy_error_log(e)


    return oneCalc,ID


###############################################################################################################

def _oneRun(__submitNodeName__,oneCalc,ID,execPrefix='',execPostfix='',engine='espresso',calcType=None,executable=None,execPath=None,nextCalc=None,nextConf=None,exit_on_error=True):
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

    logging.debug('entering _oneRun')

    os.chdir(os.path.abspath(oneCalc['_AFLOWPI_FOLDER_']))
    AFLOWpi.prep._forceGlobalConfigFile(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'../AFLOWpi/CONFIG.config'))

    clusterType = AFLOWpi.prep._ConfigSectionMap("cluster","type")

    try:
        '''save the status of the calcs'''
        if '__status__' not in list(oneCalc.keys()):
            oneCalc['__status__']=collections.OrderedDict()

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
    # '''set the completed status to false in the case of a looping set of calcs'''
    home = os.path.abspath(os.curdir)
    rd ='./'

    '''do the setup'''

    oneCalc,ID = AFLOWpi.run._qe__pre_run(oneCalc,ID,calcType,__submitNodeName__,engine)


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
            executable = AFLOWpi.run._getExecutable(engine,calcType)

            executable = os.path.join('./',executable)
        else:
            execPath = execPath
            executable=''

        if os.path.exists(os.path.join(rd,engine)):
            command='%s %s%s %s < %s  >  %s' % (execPrefix,execPath,executable,execPostfix,ri,ro)
            logging.error(command)
        else:
            command='%s %s%s %s < %s  >  %s' % (execPrefix,execPath,executable,execPostfix,ri,ro)

    except Exception as e:
        _fancy_error_log(e)
###############################################################################################################      
    try:	
        logging.info("starting %s in %s" % (command, home))
        print(('in %s:\n'%home))
        print(("starting %s\n"%command))

        if clusterType.upper() in ['PBS','SLURM','UGE']:
            job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True,cwd=rd)
            a = time.time()

            job.communicate()               

#################################################################################################################     
#################################################################################################################    
#################################################################################################################    
###CONSIDER MOVING OUTSIDE
#################################################################################################################    
#################################################################################################################    
#################################################################################################################     
            if clusterType.upper() in ['PBS','SLURM','UGE']:
                AFLOWpi.run._restartPW(oneCalc,ID,a,__submitNodeName__) 

            if job.returncode==0:
                try:
                    if calcType in ['scf']:
                        pass
                    elif calcType in ['emr','nmr','hyperfine','gvectors']:
                        __restartGIPAW(oneCalc,ID,a,__submitNodeName__)
                except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)

            if job.returncode==1 or job.returncode==255:
                logging.warning('%s in %s returned exit code: %s'%(ID,oneCalc['_AFLOWPI_FOLDER_'],job.returncode))

                try:
                    if oneCalc['__status__']['Error']!='None':
                        oneCalc['__status__']['Error']+='Exited Status: %s' % 1
                    else:
                        oneCalc['__status__']['Error']='Exit Status: %s' % 1
                except:
                    oneCalc['__status__']['Error']='Exited Status %s' % 1 
                AFLOWpi.prep._saveOneCalc(oneCalc,ID)                                
                if exit_on_error==True:
                    err_msg='Encountered Error in execution of %s during step #%02d Exiting'%(executable,oneCalc['__chain_index__'])
                    logging.error(err_msg)
                    AFLOWpi.prep._from_local_scratch(oneCalc,ID)
                    sys.exit(0)

            elif job.returncode==2:
                oneCalc['__status__']['Error']='FAILED TO CONVERGE %s' % job.returncode
                AFLOWpi.prep._saveOneCalc(oneCalc,ID)
                if exit_on_error==True:
                    err_msg='Encountered Error in execution of %s during step #%02d Exiting'%(executable,oneCalc['__chain_index__'])
                    logging.error(err_msg)
                    AFLOWpi.prep._from_local_scratch(oneCalc,ID)
                    sys.exit(0)

            elif job.returncode > 2 and job.returncode <= 127:				
                #empty text file is created in the AFLOWpi folder
                #is created and the keys of the failed to converge
                #calculations is written in the file
                oneCalc['__status__']['Error']='QE ERR %s' % job.returncode
                AFLOWpi.prep._saveOneCalc(oneCalc,ID)
                if exit_on_error==True:
                    err_msg='Encountered Error in execution of %s during step #%02d Exiting'%(executable,oneCalc['__chain_index__'])
                    logging.error(err_msg)
                    AFLOWpi.prep._from_local_scratch(oneCalc,ID)
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
#################################################################################################################     
#################################################################################################################     
#################################################################################################################     
###CONSIDER MOVING OUTSIDE
#################################################################################################################     
#################################################################################################################     
#################################################################################################################      
        logging.info("finished %s in %s" % (command, rd))
        print(("finished %s\n" % (command)))

        if AFLOWpi.prep._ConfigSectionMap('prep','save_dir') != '':
            try:

                AFLOWpi.retr._moveToSavedir(os.path.join(rd,ro))
                AFLOWpi.retr._moveToSavedir(os.path.join(rd,ri))
                logs = glob.glob('%s/*%s*.log' % (rd,oneCalc['_AFLOWPI_PREFIX_'][1:]))
                for log in logs:
                    AFLOWpi.retr._moveToSavedir(log)
                pdf = glob.glob('%s/*%s*.pdf' % (rd,oneCalc['_AFLOWPI_PREFIX_'][1:]))
            except Exception as e:
                _fancy_error_log(e)

        logging.debug('STARTING DATA')
        if oneCalc['_AFLOWPI_PREFIX_'][1:]==ID:
            logging.debug('RECORDING DATA')
            try:

                AFLOWpi.retr.getCellOutput(oneCalc,ID)

            except Exception as e:
                pass
###############################################################################################################
    except KeyboardInterrupt:
        print('Got Keyboard Exit signal, exiting')
        logging.debug('Got Keyboard Exit signal, exiting')
        sys.exit(0)
    except OSError as e:
            print(("ERROR! "+str(e)))
            logging.error("ERROR! "+str(e))
    except ValueError as e:
            print(("ERROR! "+str(e)))
            logging.error("ERROR! "+str(e))
    except subprocess.CalledProcessError as e:
            if e.returncode==1:
                print("exit code==1")

#################################################################################################################     
#################################################################################################################     
#################################################################################################################     
###CONSIDER MOVING OUTSIDE
#################################################################################################################     
#################################################################################################################     
#################################################################################################################      
    os.chdir(home)		




    logging.debug('exiting _oneRun')





def _onePrep(oneCalc,ID,execPrefix=None,execPostfix=None,engine='espresso',calcType='',alt_ID=None,exit_on_error=True):
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

    logging.debug('entering _onePrep')

    if execPrefix == None or execPrefix ==  '':
            execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
    else:
        execPrefix=''
    if execPostfix == None or execPostfix ==  '':

        execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
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
        AFLOWpi.run._oneRun''' % oneCalc['__execCounterBkgrd__']

        execString+='''(__submitNodeName__,oneCalc,%s,execPrefix='%s',execPostfix='%s',engine='%s',calcType='%s',exit_on_error=%s)
        oneCalc['__execCounter__']+=1
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
        '''  % (write_ID,execPrefix,execPostfix,engine,calcType,exit_on_error)

        oneCalc['__execCounterBkgrd__']+=1

        execFileString = os.path.abspath(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_'+ID+'.py'))
        logging.debug('ID: %s' % ID)


        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',execString)
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"'%s'" % ID)


    except Exception as e:
        print(e)
        AFLOWpi.run._fancy_error_log(e)


    logging.debug('exiting _onePrep')
################################################################################################################

################################################################################################################

def _bands_pp(__submitNodeName__,oneCalc,ID):
    
    nspin = int(AFLOWpi.scfuj.chkSpinCalc(oneCalc,ID))
    execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

    if nspin==2:
        up = oneCalc,'','bands_up',ID
        dn = oneCalc,'','bands_dn',ID

        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='bands_up',executable='bands.x',execPath='./bands.x')		


        AFLOWpi.run._PW_bands_fix(oneCalc,ID+'_up')

        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='bands_dn',executable='bands.x',execPath='./bands.x')		


        AFLOWpi.run._PW_bands_fix(oneCalc,ID+'_dn')

    else:
        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='bands',executable='bands.x',execPath='./bands.x')		

        AFLOWpi.run._PW_bands_fix(oneCalc,ID)
    


def _makeInput(oneCalc,engine,calcType,ID=''):
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
        prefix = AFLOWpi.retr._prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])
    except:
        prefix = oneCalc['_AFLOWPI_PREFIX_']

    temp_dir= AFLOWpi.prep._get_tempdir()
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
        except Exception as e:
            nbnd='30'

        stringDict={'dos':
                        """ &dos
       prefix='%s'
       outdir='%s'
       DeltaE=0.05
       degauss=0.01
       fildos='%s_dos.dat'

      /
""" % (prefix,temp_dir,ID),
                            'pdos':"""  &projwfc
       prefix='%s'
       DeltaE=0.05
!       Emax=20
!       Emin=-20
       outdir='%s'
       degauss=0.001
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

                    'bands_dn':""" &bands
     prefix='%s'
     outdir='%s'
     spin_component=2
     filband='./%s_dn_band_data.out'
 /
 """ %  (prefix,temp_dir,ID),

                    'bands_up':""" &bands
     prefix='%s'
     outdir='%s'
     spin_component=1
     filband='./%s_up_band_data.out'
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
 
        
            
################################################################################################################

################################################################################################################

# def _fixconverge(ID,folder):
#     logging.debug('entering __fixconverge')
#     logging.warning('problem converging %s/%s.in' % (folder,ID))
#     logging.debug('exiting __fixconverge') 


################################################################################################################


# def _resubmit(ID,oneCalc,__submitNodeName__,queue,calcName):

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

#         except Exception as e:
#             AFLOWpi.run._fancy_error_log(e)
#             __resubmit(IDsplit,oneCalc,__submitNodeName__,queue,calcName)


################################################################################################################



# def _lower_mixing(oneCalc,ID):


#     inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
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
#             inputString=AFLOWpi.retr._splitInput(inputDict)
#             oneCalc['_AFLOWPI_INPUT_']=inputString

#             outFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in')
#             with open(outFileName,'w') as outFile:
#                 outFile.write(inputString)
#     except:
#         pass
#     AFLOWpi.prep._saveOneCalc(oneCalc,ID)
    
#     return oneCalc,ID
    

import AFLOWpi
import os
import glob
import shutil
import re
import numpy 


def _grab_elastic_generated_inputs(oneCalc,ID):
    glob_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'Structures_ESPRESSO')
    input_files = glob.glob(glob_path+'/*.in')

    return input_files

def _prep_elastic(oneCalc,ID,eta_max=0.005,num_dist=49,use_stress = True,order=2):

    if ID in oneCalc['prev']:
        return

    try:
        os.rm(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'ElaStic_PW.in'))
    except:
        pass

    split_input = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    pos = AFLOWpi.retr.detachPosFlags(split_input['ATOMIC_POSITIONS']['__content__'])[0]
    split_input['ATOMIC_POSITIONS']['__content__']=pos

    oneCalc['_AFLOWPI_INPUT_'] = AFLOWpi.retr._joinInput(split_input)
    
    write_in = oneCalc['_AFLOWPI_INPUT_']
    elastic_qe_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'ElaStic_PW.in')
    with open(elastic_qe_in,'w') as eqfo:
        eqfo.write(write_in)

    elastic_ID = ID+'_ElaStic'

    if use_stress:
        calc_meth=2
    else:
        calc_meth=1

    elastic_in_string='''%s
%s
%s
%s'''%(calc_meth,order,eta_max,num_dist)

    elastic_in_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%elastic_ID)
    with open(elastic_in_file,'w') as eifo:
        eifo.write(elastic_in_string)

    os.system('ElaStic_Setup_ESPRESSO<%s'%elastic_in_file)



def _copy_qe_out_file(oneCalc,ID):
    stripped_step_ID='_'.join(ID.split('_')[:-1])
    dist_num = ID.split('_')[0]

    orig_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out'%ID)
    new_path  = os.path.join('../../%s/%s/'%(dist_num,stripped_step_ID),'%s.out'%(stripped_step_ID))

    try:
        shutil.copy(orig_path,new_path)
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

def _grab_el_plot_data(file_path):
    fit_data_by_order={}

    with open(file_path,'r') as dfo:
        dfs = dfo.read()

    per_fit_order = re.findall('#.*(\d+).*\n((?:[-0-9.\s]*\n)*)',dfs)
    for order,fit_data in per_fit_order:
        fit_array=[]
        for i in fit_data.split('\n'):
            to_float = list(map(float,i.split()))
            if len(to_float)!=0:
                fit_array.append(to_float)

        fit_data_by_order[order]=fit_array

    return fit_data_by_order

def _find_flat_region(data):
    eta = []
    val = []
    #Sort the data into two arrays. one for eta and one for the values
    for k,v in data:
        eta.append(k)
        val.append(v)

    #set the max variance from one step to the next for the plateau region.
    #if the variance from one step to the next is greater than the value
    #below then that is considered the end of the plateau.
    variance=10.0

    plateau_index=0
    #reverse the order or eta and values to start with lowest eta first
    val = [x for x in reversed(val)]
    eta = [x for x in reversed(eta)]

    #find difference between one point and the next
    grad = numpy.gradient(val).tolist()
    found=False
    for x in range(len(grad)):
        #if we find a values between one eta and the next is less than
        #our designated variance threshold then start looking for the end
        #if it. If we haven't found a plateau and the variance is greater
        #than the threshold keep looking for the plateau. If we have found
        #the plateau and we find the end of it then record the index of the
        #end of the plateau and use the eta at that index as our max eta
        if grad[x]<variance:
            if found==False:
                found=True

            plateau_index=x
        else:
            if found==False:
                continue
            else:
                break

    print(('plateau ends at index:',plateau_index))

    return eta[plateau_index]


def _pp_elastic(oneCalc,ID,order=2,use_stress=True):

    res_analyze=''
    result_input_file = ''
    result_output_file = ''
    if use_stress:
        os.system('ElaStic_Analyze_Stress')
        if order==2:
            res_analyze = 'ElaStic_Result_Stress_2nd'
            result_input_file = 'ElaStic_2nd.in'
            result_output_file = 'ElaStic_2nd.out'
        if order==3:
            res_analyze = 'ElaStic_Result_Stress_3rd'
            result_input_file = 'ElaStic_3rd.in'
            result_output_file = 'ElaStic_3rd.out'
    else:
        os.system('ElaStic_Analyze_Energy')

        if order==2:
            res_analyze = 'ElaStic_Result_Energy_2nd'
            result_input_file = 'ElaStic_2nd.in'
            result_output_file = 'ElaStic_2nd.out'
        if order==3:
            res_analyze = 'ElaStic_Result_Energy_3rd'
            result_input_file = 'ElaStic_3rd.in'
            result_output_file = 'ElaStic_3rd.out'



    svs_folder = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'Stress-vs-Strain')

    fit_files = glob.glob(svs_folder+'/*_d1S.dat')


    files_sep_dist={}
    results_input_str=''
    for data_file in fit_files:
        distortion_ID = os.path.basename(data_file).split('_')[0]
        try:
            files_sep_dist[distortion_ID].append(data_file)
        except:
            files_sep_dist[distortion_ID]=[data_file]
    
    
    for distortion_ID,fit_files in list(files_sep_dist.items()):
        order_list=[]
        m_eta_list=[]
        for data_file in fit_files:
            distortion_ID = os.path.basename(data_file).split('_')[0]
            fit_data = AFLOWpi.run._grab_el_plot_data(data_file)
            #set a global max eta for all fittings
            global_max_eta=0.0
            #set the order of the polynomial to fit to.
            poly_fit_choice = 0

            if use_stress:
                fit_index=1
            else:
                fit_index=2
            #iterate through all orders of poly fit data
            for poly_order,data in list(fit_data.items()):
                #find max eta for one order of poly fit
                max_eta = AFLOWpi.run._find_flat_region(data)
                #check if the max eta is greater or equal to the current global
                #max eta. If the value is equal then prefer a higher order
                #polynomial fit to a lower one.
                if max_eta>=global_max_eta:
                    global_max_eta=max_eta
                    fit_index=poly_order
            #append the values for max eta and poly fit order 
            #to a list to later be used to create the the ElaStic_2nd.in
            #file.
            order_list.append(fit_index)
            m_eta_list.append(global_max_eta)
        #write the entry in the file
        results_input_str+=distortion_ID+' '+' '.join(map(str,m_eta_list))+'\n'
        results_input_str+=' '*len(distortion_ID)+' '.join(map(str,order_list))+'\n'

    #write the input for ElaStic_Result_Stress_2nd
    res_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],result_input_file)
    with open(res_input,'w') as rifo:
        rifo.write(results_input_str)

    #run ElaStic_Result_Stress_2nd or whatever result analyzer code
    #that is to be run.
    os.system(res_analyze)
    #move results to file with unique name
    res_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],result_output_file)
    new_res_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_elastic.out'%ID)

    shutil.move(res_output,new_res_output)

import AFLOWpi
import os
import re
import copy 

def _collect_fd_field_forces(oneCalc,ID,for_type='raman'):
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
        num=list(range(19))
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

def _pull_polarization(oneCalc,ID):
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



    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        return
    e_pol_split = [float(x.split()[-1]) for x in ele_block.split('\n') if len(x.strip())!=0]
    e_ion_split = [float(x.split()[-1]) for x in ion_block.split('\n') if len(x.strip())!=0]
    
    
    e_tot_split = [e_pol_split[0],e_pol_split[1],e_pol_split[2],]

    epol_out_string=''

    epol_out_string+='%25.15e%25.15e%25.15e\n'%(float(e_tot_split[0]),float(e_tot_split[1]),float(e_tot_split[2]),)

    force_postfix='.'.join(ID.split('.')[1:])
    force_postfix=force_postfix.split('_')[0]

    with open(os.path.join(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']),'epol.%s'%force_postfix),'w+') as out_file_obj:
        out_file_obj.write(epol_out_string)


def _pull_eps_out(oneCalc,ID):
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
    eps_regex=re.compile('Dielectric tensor\n.*\n\s*([0-9\s.-]*)\n')
    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_epol.out'%ID),'r') as out_file_obj:
            out_string = out_file_obj.read()

        eps_string=eps_regex.findall(out_string)[-1]
        
        
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    return eps_string

def _pull_born_out(oneCalc,ID):
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
            born_string+='         %s\n'%num_index
            born_string+=born_charges[i]

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    return born_string

def _gen_fd_input(oneCalc,ID,for_type='raman',de=0.003):
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

    forces = AFLOWpi.run._collect_fd_field_forces(oneCalc,ID,for_type=for_type)
    header+=forces
    
    finite_field_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_%s.in'%(ID,for_type))
    with open(finite_field_input,'w') as force_file:
        force_file.write(header)

def _swap_scf_inputs(oneCalc,ID):
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

def _setup_raman(calc_subset,orig_oneCalc,orig_ID,field_strength=0.001,field_cycles=3,for_type='raman'):
    """
    Sets up the input files and the command to run the finite field calculations.
    
    Arguments:
          calc_subset (dict): dictionary of dictionaries of the FD_PHONON calculations 
          orig_oneCalc (dict): oneCalc of one calculation the main calc set

    Keyword Arguments:
         feild_strength (float): applied electric field strength
         field_cycles (int) number of berry phase cycles
         for_type (str): which type of FD input do we want to generate files for. 'raman','born','eps'

    Returns:
         None
          
    """

    input_string=orig_oneCalc['_AFLOWPI_INPUT_']
    subset_keys = list(calc_subset.keys())

    in_dict=AFLOWpi.retr._splitInput(input_string)

    #if possibly a metal don't do LOTO

    if 'degauss' not in  list(in_dict['&system'].keys()):
        pass
    else:
        bandgap_type = AFLOWpi.retr.gap_type(orig_oneCalc,orig_ID)
        if bandgap_type not in ['p-type','n-type','insulator']:
            return
        else:
            if 'degauss' in list(in_dict['&system'].keys()):
                del in_dict['&system']['degauss']
            else:
                pass

            if "occupations" in list(in_dict['&system'].keys()):
                val=eval(in_dict['&system']['occupations'].lower())
                if val=="smearing":
                    del in_dict['&system']['occupations']
                else:
                    pass
            else:
                pass

        input_string=AFLOWpi.retr._joinInput(in_dict)


    #copy for when we append finite field calcs to the calclog
    calc_subset_extended=copy.deepcopy(calc_subset)

    if for_type=='raman':
        val_type='force'
        num=list(range(19))
    elif for_type=='born':
        val_type='force'
        num=[0,2,4,6,1,3,5]
    elif for_type=='epol':
        val_type='epol'
        num=[0,2,4,6,1,3,5]

    num_calcs=19
    subset_size=len(calc_subset)
    
    for index in range(len(num)):        
        fd_input = AFLOWpi.run._field_factory(input_string,field_strength=field_strength,for_type=for_type,field_cycles=field_cycles,dir_index=num[index])

        subset_key_index=num[index]%subset_size
        exec_counter_index = (index/subset_size)+1
        ID = subset_keys[subset_key_index]
        oneCalc=calc_subset[ID]

        #set prefix in finite field input to that of the FD phonon calc's input
        split_input = AFLOWpi.retr._splitInput(fd_input)
        split_input['&control']['prefix']='"%s"'%calc_subset[ID]['_AFLOWPI_PREFIX_']
        fd_input = AFLOWpi.retr._joinInput(split_input)

        FD_ID='finite_field.%02d'%num[index]
        file_path_fd_field_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%FD_ID)
    

        with open(file_path_fd_field_in,'w') as fd_in_obj:
            fd_in_obj.write(fd_input)


        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"""oneCalc_ff = AFLOWpi.run._swap_scf_inputs(oneCalc,'%s')""" % FD_ID)       
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','''AFLOWpi.prep._to_local_scratch(oneCalc_ff,'%s')'''%FD_ID)


        execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
        run_command="""if oneCalc['__execCounter__']<=%s:
        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s',execPrefix='%s',execPostfix='-northo 1',engine='espresso',calcType='scf',exit_on_error=True)
        oneCalc['__execCounter__']+=1
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
"""%(str(exec_counter_index),FD_ID,execPrefix)

        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',run_command)       

#        AFLOWpi.run._onePrep(oneCalc,ID,engine='espresso',calcType='scf',execPrefix=execPrefix,execPostfix=' -northo 1 ',alt_ID=FD_ID,)


        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"""AFLOWpi.run._pull_forces(oneCalc_ff,'%s')""" % FD_ID)       
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN',"""AFLOWpi.run._pull_polarization(oneCalc_ff,'%s')""" % FD_ID)       


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

#    AFLOWpi.prep.updatelogs(calc_subset_extended,chain_logname,runlocal=True)

def _field_factory(input_str,field_strength=0.001,for_type='raman',field_cycles=3,dir_index=0):
    """
    Modifies the QE pwscf input for the finite field calc of a given index
    
    Arguments:
         input_str (str): string of QE pwscf input file

    Keyword Arguments:
         feild_strength (float): applied electric field strength
         field_cycles (int) number of berry phase cycles
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


    in_dict=AFLOWpi.retr._splitInput(input_str)
    in_dict['&control']['lelfield']='.TRUE.'
    in_dict['&control']['calculation']="'scf'"
    in_dict['&control']['nberrycyc']=field_cycles
    in_dict['&control']['tprnfor']='.TRUE.'
    in_dict['&control']['prefix']="'_finite_field.%02d'"%dir_index

    input_list=[]

#    for i in range(len(field_dir)):
    for j in range(len(field_dir[dir_index])):
        in_dict['&electrons']['efield_cart(%s)'%(j+1)]=field_dir[dir_index][j]*field_strength

    new_input = AFLOWpi.retr._joinInput(in_dict)

        
    return new_input

import AFLOWpi
import numpy 
import scipy
import os
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot


def _gen_smear_conv_calcs(oneCalc,ID,num_points=4,smear_type='mp',smear_variance=0.3,calc_type='scf'):

    input_dict = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])

    try:
        initial_degauss = input_dict['&system']['degauss']
    except:
        initial_degauss = '0.01'

        input_dict['&system']['occupations']='"smearing"'
        if 'smearing' not in list(input_dict['&system'].keys()):
            input_dict['&system']['smearing']='"%s"'%smear_type

    input_dict['&system']['degauss']=initial_degauss

    fl_degauss = float(initial_degauss)
    ub_degauss = fl_degauss*(1.0+smear_variance)
    lb_degauss = fl_degauss*(1.0-smear_variance)
    
    ub_degauss = fl_degauss*(5.0)
    lb_degauss = fl_degauss*(0.5)

    ub_degauss = fl_degauss+0.005
    lb_degauss = fl_degauss-0.005
    degauss_vals = numpy.linspace(lb_degauss,ub_degauss,num_points)


    input_files = []

    for val in degauss_vals:
        input_dict['&system']['degauss'] = val
        input_dict['&control']['calculation'] = "'%s'"%calc_type
        one_input = AFLOWpi.retr._joinInput(input_dict)

        input_files.append(one_input)


    return input_files



def _extrapolate_smearing(oneCalc,ID):

    folder = oneCalc['_AFLOWPI_FOLDER_']
    calc_log = os.path.join(folder,'SMEARING','AFLOWpi','calclogs','step_01.log')

    calcs=AFLOWpi.prep._load_log_from_filename(calc_log)

    de_set=[]
    smear_vals=[]

    calcs = AFLOWpi.retr.grabEnergyOut(calcs)
    for k,v in list(calcs.items()):

        oc_in_dict = AFLOWpi.retr._splitInput(v['_AFLOWPI_INPUT_'])
        smear_type = oc_in_dict['&system']['smearing']
        smear_val = float(oc_in_dict['&system']['degauss'])
        en         = float(v['Energy'])
#        force=AFLOWpi.retr.getForce(v,k,string=False)
        if en != 0.0:
#            de_set.append(force)
            de_set.append(en)
            smear_vals.append(smear_val)


    de_set = numpy.asarray(de_set)
    #CHECK THIS
    smear_vals = numpy.asarray(smear_vals)
#    (x_0,x_1) = 
#    de_set=numpy.abs(de_set)
    max_de = numpy.amax(smear_vals)
    min_de = numpy.amin(smear_vals)

    if min_de<=0:
        min_end=2.0*min_de
    else:
        min_end=-2.0*min_de

    if max_de<=0:
        max_end=-2.0*max_de
    else:
        max_end= max_de*2.0

    extrap_smear= numpy.linspace(0.0,max_end,200)
#    interp_en = interp_en.tolist()
#    interp_en.append(0.0)
#    interp_en = numpy.asarray(interp_en)

    params=[smear_vals,de_set,]

#    extrap_smear = numpy.polyval(numpy.polyfit(en_set,smear_vals,2),interp_en)

#    en_data_diff=numpy.abs(numpy.gradient(params[1]))
#    en_data_diff=numpy.gradient(params[1])
    fit = numpy.polyfit(params[0],params[1],2)
    fit_en =numpy.polyval(fit,extrap_smear)
#    en_deriv = numpy.polyfit(params[0],en_data_diff,1)
#    print fit
    en_deriv=numpy.polyder(fit)
#    en_deriv=fit

    print(en_deriv)
#    en_deriv=fit.deriv()
#    fir1d= numpy.polyfit(params[0],params[1],1)
    interp_en =numpy.polyval(en_deriv,extrap_smear)

#    interp_en = en_deriv(extrap_smear)
#    at_zero = en_deriv[0]
    at_zero=en_deriv[0]/en_deriv[1]
    print(at_zero)
    print(de_set)
#    print en_data_diff
    inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

    if at_zero<0.001:
        #if we don't need smearing get rid of it.
        try:
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','degauss',at_zero,del_value=True)
        except:
            pass
        try:
            #CHECK THIS!!!!
            smu = eval(inputDict['&system']['occupations'])
            if smu == 'smearing':
                oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','occupations','"smearing"',del_value=True)

        except:
            pass
        try:
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','smearing',smear_type,del_value=True)
        except: 
            pass
        

    else:
        #if we need smearing
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','degauss',at_zero,del_value=False)
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','occupations','"smearing"',del_value=False)
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','smearing',smear_type,del_value=False)


    return_input = AFLOWpi.retr._joinInput(inputDict)
    oneCalc['_AFLOWPI_INPUT_']=return_input


    matplotlib.rc("font", family="serif")      #to set the font type
    matplotlib.rc("font", size=10)             #to set the font size

    matplotlib.pyplot.plot(smear_vals,de_set,"b.")
    matplotlib.pyplot.plot(extrap_smear,fit_en,"r")
#    matplotlib.pyplot.plot(smear_vals,en_data_diff,"b.")
#    matplotlib.pyplot.plot(extrap_smear,interp_en,"r")
    matplotlib.pyplot.ylabel("|$\Delta$E|")
    matplotlib.pyplot.xlabel("degauss (Ry)")
    subdir=oneCalc['_AFLOWPI_FOLDER_']
    fileplot = os.path.join(subdir,'SMEARING_%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID))

    matplotlib.pyplot.savefig(fileplot)

    return oneCalc,ID



import AFLOWpi
import os
import glob
import subprocess
import collections
import logging
import sys
import re
import time 

def _get_cluster_submit_command():
    '''set parent processID to 0 so they'll be different when the script starts'''
    clusterType = AFLOWpi.prep._ConfigSectionMap("cluster","type").upper()
    submitCommand='qsub'
    if clusterType.upper()=='SLURM':
        submitCommand='sbatch'  
    elif clusterType.upper()=='UGE':
            submitCommand='qsub'


    return submitCommand
                    
def _check_and_submit(submit_file,sub_command):
    if os.path.exists(submit_file):
        qf=re.sub('.submit','.qsub',submit_file)
        os.system('%s %s'%(sub_command,qf))
        directory=os.path.dirname(submit_file)
        ID=os.path.basename(submit_file)[1:-7]
        try:
            with open("./submitted.log", "a") as submitted_log_file:
                submitted_log_file.write('%s\n'%ID)
        except:
            with open("./submitted.log", "w") as submitted_log_file:
                submitted_log_file.write('%s\n'%ID)

        logging.info('Daemon submitted %s in %s'%(ID,directory))
        os.remove(submit_file)

def _generate_submission_daemon_script(calcs):
    workdir=AFLOWpi.prep._ConfigSectionMap('prep','work_dir')
    for ID,oneCalc in list(calcs.items()):
        project=oneCalc['PROJECT']
        calc_set=oneCalc['SET']
        break
    configFile = os.path.join(workdir,project,calc_set,'AFLOWpi','CONFIG.config')
    submit_daemon_dir = os.path.join(workdir,project,calc_set,'AFLOWpi','submission_daemon')
    if not os.path.exists(submit_daemon_dir):
        os.mkdir(submit_daemon_dir)

    daemon_file_name = os.path.join(submit_daemon_dir,'submit_daemon.py')

    daemon_file_string='''
import AFLOWpi
import logging
import time

logging.basicConfig(filename='../LOG.log',format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG)
configFile='''

    daemon_file_string+=repr(configFile)
    daemon_file_string+='''

AFLOWpi.prep._forceGlobalConfigFile(configFile)
AFLOWpi.run._run_submission_check()
AFLOWpi.run._restart_submission_daemon('submit_daemon.py')'''

    print(('generating daemon script in %s'%daemon_file_name))
    with open(daemon_file_name,'w') as daemon_file_object:
        daemon_file_object.write(daemon_file_string)


    AFLOWpi.run._start_submission_daemon(daemon_file_name)

def _start_submission_daemon(daemon_file_name):
    print(('starting daemon script: %s'%daemon_file_name))    
    cur_dir=os.curdir
    os.chdir(os.path.dirname(daemon_file_name))

    subprocess.Popen('nohup python ./submit_daemon.py  &',shell=True)
    os.chdir(cur_dir)
    print('daemon script started')

def _restart_submission_daemon(file_name):

    subprocess.Popen('nohup python ./%s &'%file_name,shell=True)





def _submit_log_append(addition_list,logname='./log_list.log'):
    daemon_submission_dir = os.path.dirname(logname)
    if not os.path.exists(daemon_submission_dir):
        os.mkdir(daemon_submission_dir)
    if os.path.exists(logname):
        with open(logname,'r') as log_list_obj:
            log_list_string = log_list_obj.read()

        logs=log_list_string.split('\n')    
    else:
        logs=[]

    logs.extend(addition_list)
    logs=list(set(logs))
    log_list_string='\n'.join(logs)

    with open(logname,'w') as log_list_obj:
        log_list_obj.write(log_list_string)

def _run_submission_check():
    calc_list = AFLOWpi.run._load_submit_log()
    
    need_submitting=AFLOWpi.run._check_statuses(calc_list)

    sub_command = AFLOWpi.run._get_cluster_submit_command()
    for submit_file in need_submitting:
        AFLOWpi.run._check_and_submit(submit_file,sub_command)


def _load_submit_log():
    with open('./log_list.log','r') as log_list_obj:
        log_list_string = log_list_obj.read()
    
    logs=log_list_string.split('\n')

    calc_list=[]
    for filename in reversed(logs):
        calc_step = AFLOWpi.run._get_potential_sub_locs(filename)
        calc_list.append(calc_step)

    return calc_list

def _get_potential_sub_locs(filename):
    with open(filename,'r') as log_file_obj:
        log_file_str = log_file_obj.read()

    oneCalc_locs = [loc for loc in log_file_str.split('\n') if len(loc.strip())]
    calcs={}
    for oneCalc_file in oneCalc_locs:
        ID= os.path.basename(oneCalc_file)[1:-8]
        folder = os.path.dirname(oneCalc_file)
        if os.path.exists(oneCalc_file):
            oneCalc=AFLOWpi.prep._loadOneCalc(folder,ID)
            calcs[ID]=oneCalc
        else:
            calcs[ID]={'_AFLOWPI_FOLDER_':folder}
            
    return calcs


def _check_statuses(calc_list):

    keep_alive=False
    chain_status={}
    submission_check_list=[]
    calc_list_by_chain = AFLOWpi.run._sort_by_chain(calc_list)

    for prefix in list(calc_list_by_chain.keys()):
        chain = calc_list_by_chain[prefix]

        for ID in list(chain.keys()):
            try:
                completed=chain[ID]['__status__']['Complete']
                error=chain[ID]['__status__']['Error']
                if completed!=True:
                    if error!='None':
                        break
                    else:
                        keep_alive=True
                        dir_name = chain[ID]['_AFLOWPI_FOLDER_']
                        check_submit_file=os.path.join(dir_name,'_%s.submit'%ID)
                        submission_check_list.append(check_submit_file)
            except:
                dir_name = chain[ID]['_AFLOWPI_FOLDER_']
                check_submit_file=os.path.join(dir_name,'_%s.submit'%ID)
                submission_check_list.append(check_submit_file)
                keep_alive=True

    '''check the calc_list log again just in case it got updated while the daemon was sleeping'''
    time.sleep(60)
    calc_list_new = AFLOWpi.run._load_submit_log()
    if calc_list_new!=calc_list:
        keep_alive=True

    if keep_alive==False:
        logging.info('All calculations in set are finished. Stopping daemon')
        print('All calculations in set are finished. Stopping daemon')
        sys.exit(0)
    else:
        return submission_check_list


def _sort_by_chain(calc_list):
    by_chain={}

    for calc_set in calc_list:
        for ID,oneCalc in list(calc_set.items()):        
            try:
                prefix=oneCalc['_AFLOWPI_PREFIX_']
            except:
                prefix='_'+ID.split('_')[0]+'_01'
            try:
                by_chain[prefix][ID]=oneCalc
            except:
                by_chain[prefix]=collections.OrderedDict()
                by_chain[prefix][ID]=oneCalc

    return by_chain

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
    except Exception as e:
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


            except Exception as e:
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
      except Exception as e:
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
        except Exception as e:
    	    AFLOWpi.run._fancy_error_log(e)
    #remove the fd.x output files in case we do another phonon calc
    try:
	    globpath=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'FD_PHONON')
	    phil = glob.glob(globpath+'/displaced*.in')
	    for i in phil:
		    os.remove(i)
    except Exception as e:
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

    except Exception as e:
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
        qs=[list(map(float,k.split()[1:-1])) for k in eig]
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

    kpts,data,summed_weights,species = AFLOWpi.run._get_ph_weights(appdos_fn)

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

    data = numpy.asarray([list(map(float,x.split())) for x in data[1:] if len(x.strip())!=0])
    
    kpts=data[:3]
    data=data[:,3:]

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

    return kpts,data,summed_weights,species

