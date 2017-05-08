# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
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
import datetime
import cPickle
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

    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.prep._addToBlock(oneCalc,ID,'CLEANUP','\nimport os\nos.system("rm -rf %s/_*")\n' % oneCalc['_AFLOWPI_FOLDER_'])            
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

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

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
            print '%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%'
            try:
                print e
                for errorMSG in errorList:
                    print errorMSG

            except Exception,e:
                print e

            print '%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%DEBUG%%%'
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
    except Exception,e:
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

    print 'EMR NOT IMPLEMENTED'
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
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)

    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='emr')
        except Exception,e:
            _fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='emr')
        except Exception,e:
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
    print 'HYPERFINE NOT IMPLEMENTED'
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
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)


    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='hyperfine')
        except Exception,e:
            _fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='hyperfine')
        except Exception,e:
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

    print 'NMR NOT IMPLEMENTED'
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
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)



    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='nmr')
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='nmr')
        except Exception,e:
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

    print 'GVECTORS NOT IMPLEMENTED'
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
        print 'ERROR: gipaw not found in %s please check your config file EXITING' % gipawdir
        raise SystemExit
    AFLOWpi.prep.totree(gipawExLoc, calcs,rename=None,symlink=symlink)


    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='gvectors')
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='gvectors')
        except Exception,e:
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
    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',exit_on_error=exit_on_error)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf')
        except Exception,e:
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
    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='dos')
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='dos')
        except Exception,e:
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
    for ID,oneCalc in calcs.iteritems():
        try:
            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='pdos')
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
        try:
            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='pdos')
        except Exception,e:
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
    for ID,oneCalc in calcs.iteritems():
#        try:
#            AFLOWpi.run._onePrep(oneCalc,ID,execPrefix=execPrefix,execPostfix=' ',engine='espresso',calcType='bands')
#        except Exception,e:
#            AFLOWpi.run._fancy_error_log(e)
#        try:
#            AFLOWpi.run._testOne(ID,oneCalc,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='bands')
#        except Exception,e:
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
            except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)

            for files in enginePath:

                if os.path.exists(os.path.join(engineDir,files))==False:
                    print '%s not found. Check your config file to make sure engine_dir path is correct and that %s is in that directory..Exiting' %(files,files)
                    logging.error('%s not found. Check your config file to make sure engine_dir path is correct and that %s is in that directory..Exiting' %(files,files))
                    raise SystemExit

                if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
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
                CONFIG_FILE='../AFLOWpi/CONFIG.config'
                calcs[ID]['_AFLOWPI_CONFIG_']=CONFIG_FILE
                calcs[ID]['calcType']=calcType                

        except Exception,e:
            _fancy_error_log(e)
        #############################################################################################################
        clusterType = AFLOWpi.prep._ConfigSectionMap("cluster",'type')
        if clusterType=='':
            clusterType=''
        try:
            if 'firstCalcList' not in globals().keys() or holdFlag==False:

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

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)

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



    except Exception,e:
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
    except Exception,e:
        _fancy_error_log(e)

    try:
        engineDir  = AFLOWpi.prep._ConfigSectionMap("prep",'engine_dir')
        if os.path.isabs(engineDir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            engineDir =  os.path.join(configFileLocation, engineDir)
    except Exception,e:
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
    except Exception,e:
        _fancy_error_log(e)

    try:
        stepsAsJobs = AFLOWpi.prep._ConfigSectionMap("cluster","steps_as_jobs").lower()
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
    except Exception,e:
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
                
            except Exception,e:
                _fancy_error_log(e)
            
            qsubFilename = os.path.abspath(os.path.join(folder,'_'+submit_ID+'.qsub'))

            calcName=AFLOWpi.run._get_qsub_name(oneCalc['_AFLOWPI_FOLDER_'])

            queue = ''
            if AFLOWpi.prep._ConfigSectionMap("cluster","queue") !='':
                queue = '-q %s ' %  AFLOWpi.prep._ConfigSectionMap("cluster","queue")
        except Exception,e:
            _fancy_error_log(e)
        try:
            
            command  = "ssh -o StrictHostKeyChecking=no %s '%s %s -N %s %s %s' " % (__submitNodeName__,submitCommand,queue,calcName,qsubFilename,additional_flags)
            job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True)
            job.communicate()               

            logging.info('submitted %s to the queue' % submit_ID)
        except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
            logging.info("Trying to submit job with no -N option")

            try:
                command = "ssh -o StrictHostKeyChecking=no %s '%s %s %s' " % (__submitNodeName__,submitCommand,qsubFilename,additional_flags)
                job =  subprocess.Popen(command,stderr = subprocess.PIPE,shell=True)
                job.communicate()               

            except Exception,e:
                logging.info("Trying to submit job without ssh and -N option")
                try:
                    baseID = ID.split('_')[0]
                    baseID = ID
                    qsubFilename = os.path.abspath(os.path.join(folder,'_'+baseID+'.qsub'))
                    os.system("ssh -o StrictHostKeyChecking=no %s '%s %s -N %s %s %s' " % (__submitNodeName__,submitCommand,queue,calcName,qsubFilename,additional_flags))
                except Exception,e:
                    logging.info("Trying to submit job without ssh and -N option")
                    try:
                        baseID = ID.split('_')[0]
                        baseID = ID
                        qsubFilename = os.path.abspath(os.path.join(folder,'_'+baseID+'.qsub'))
                        os.system("%s %s %s %s" % (submitCommand,queue,qsubFilename,additional_flags))
                    except Exception,e:
                        logging.info("Trying to submit job without ssh and -N option")
                        try:
                            os.system("%s %s %s %s" % (submitCommand,queue,qsubFilename,additional_flags))
                        except Exception,e:
                            _fancy_error_log(e)

    else:
        try:

            execFileString  = os.path.abspath(os.path.join(folder,'_'+ID+'.py'))
            
            execfile(execFileString)

        except KeyboardInterrupt:
            print 'Got Keyboard Exit signal, exiting'
            logging.debug('Got Keyboard Exit signal, exiting')
            sys.exit(0)
        except Exception,e:
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
    if '__submitNodeName__' not in globals().keys():
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

    print AFLOWpi.run._colorize_message(submit_message,level='ERROR',show_level=False)

    for ID,oneCalc in calcs.iteritems():
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

    except Exception,e:
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
    except Exception,e:
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
    except Exception,e:
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

        except Exception,e:
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
            restartInput = AFLOWpi.retr._joinInput(inputDict)
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.in" % ID),"w") as inputfile:
                inputfile.write(restartInput)
            oneCalc['_AFLOWPI_INPUT_']=restartInput
        except Exception,e:
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
    except Exception,e:
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

        if  '__GLOBAL_CLOCK__' not in globals().keys():
            '''only move from local scratch if this script is starting fresh'''

            AFLOWpi.prep._to_local_scratch(oneCalc,ID)
            oneCalc['__walltime_dict__']={}
            globals()['__GLOBAL_CLOCK__']=__CLOCK__
            walltime_dict['walltime'] = AFLOWpi.run._getWalltime(oneCalc,ID)         
            walltime_dict['start']=__CLOCK__
            AFLOWpi.run._writeWalltimeLog(oneCalc,ID,walltime_dict)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)
        try: 

            oneCalc['__walltime_dict__']={}
            walltime_dict={}
            walltime_dict['walltime'] = 40000000000.0         
            walltime_dict['start']=0.0

            AFLOWpi.run._writeWalltimeLog(oneCalc,ID,walltime_dict)
        except Exception,e:
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

    oneCalc =AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],ID)
    try:
        oneCalc['__status__']['Restart']+=1
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)
        AFLOWpi.run._submitJob(ID,oneCalc,__submitNodeName__,sajOverride=True)
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
    except Exception,e:
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
    except Exception,e:
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




    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

    try:
        if calcType!='scf' and calcType!=None and calcType!='custom':
                with open(os.path.join(rd,"%s_%s.in") % (ID,calcType),'w+') as PPInput:
                    PPInput.write(AFLOWpi.run._makeInput(oneCalc,engine,calcType,ID=ID))

    except Exception,e:
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
        if '__status__' not in oneCalc.keys():
            oneCalc['__status__']=collections.OrderedDict()

    except Exception,e:
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

    except Exception,e:
        _fancy_error_log(e)
###############################################################################################################      
    try:	
        logging.info("starting %s in %s" % (command, home))
        print 'in %s:\n'%home
        print "starting %s\n"%command

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
                except Exception,e:
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
        print "finished %s\n" % (command)

        if AFLOWpi.prep._ConfigSectionMap('prep','save_dir') != '':
            try:

                AFLOWpi.retr._moveToSavedir(os.path.join(rd,ro))
                AFLOWpi.retr._moveToSavedir(os.path.join(rd,ri))
                logs = glob.glob('%s/*%s*.log' % (rd,oneCalc['_AFLOWPI_PREFIX_'][1:]))
                for log in logs:
                    AFLOWpi.retr._moveToSavedir(log)
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


    except Exception,e:
        print e
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
        efermi=AFLOWpi.retr._getEfermi(oneCalc,ID,directID=False)
    except:
        efermi=0.0

    emax=40.0-0025
    emin=-40-0.0025

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
        except Exception,e:
            nbnd='30'

        stringDict={'dos':
                        """ &dos
       prefix='%s'
       outdir='%s'
 !      Emin=%s
 !      Emax=%s
       degauss=0.025
       fildos='%s_dos.dat'

      /
""" % (prefix,temp_dir,emin,emax,ID),
                            'pdos':"""  &projwfc
       prefix='%s'
       DeltaE=0.005
       outdir='%s'
       Emin=%s
       Emax=%s
       degauss=0.005
    !    kresolveddos=.true.
       filpdos='%s'

      / 
 """ % (prefix,temp_dir,emin,emax,ID),

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

#         except Exception,e:
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
    
