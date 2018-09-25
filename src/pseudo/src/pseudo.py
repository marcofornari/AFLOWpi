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
import copy
from collections import OrderedDict
import numpy
from scipy.optimize import curve_fit
import AFLOWpi.prep 
import AFLOWpi.run
import logging
import logging.handlers
import itertools as it
import copy 
import time
import sys
import random
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import glob
import shutil
import socket
import ast

def _passGlobalVar(varname,value):
    globals()[varname]=value


#def _shift_curoffs(oneCalc,ID):
#    return oneCalc,ID

def brute_test(calcs,ecutwfc,dual=None,sampling=None,constraint=None,thresh=0.001,initial_variance=0.05,grid_density=10,mult_jobs=False,conv_thresh=0.01,calc_type='relax'):

    for ID,oneCalc in calcs.iteritems():
        #set sampling and dual list from input if no variance provided
        if sampling==None or dual==None:
            in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
            sampling = [in_dict['K_POINTS']['__content__'], ]
        if dual==None:
            dual_from_input=4
            wfc_cut = float(in_dict['&system']['ecutwfc']) 
            try:
                rho_cut = float(in_dict['&system']['ecutrho']) 
                dual_from_input = int(rho_cut/wfc_cut)
                dual = [dual_from_input, ]                
            except:
                dual = [dual_from_input, ]         


        #list of the selected cutoff and sampling
#        oneCalc['_AFLOWPI_PPTEST_ECUTWFC_']=sorted(ecutwfc)
#        oneCalc['_AFLOWPI_PPTEST_SAMPLING_']=sorted(sampling)
#        oneCalc['_AFLOWPI_PPTEST_DUAL_']=sorted(dual)
        oneCalc['_AFLOWPI_PPTEST_PREV_']={}

        #tracker for the sampling and cutoff


        AFLOWpi.prep._addToBlock(oneCalc,ID,'PREPROCESSING','oneCalc,ID=AFLOWpi.pseudo._set_cutoffs(oneCalc,ID,ecutwfc=%s,dual=%s,sampling=%s)'%(ecutwfc,dual,sampling))
        AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','oneCalc,ID=AFLOWpi.pseudo._check_test_conv(oneCalc,ID,%s)'%thresh)

    AFLOWpi.pseudo.crawlingMinimization(calcs,faultTolerant=True,constraint=constraint,thresh=thresh,initial_variance=initial_variance,grid_density=grid_density,mult_jobs=mult_jobs,final_minimization=None)

def _set_cutoffs(oneCalc,ID,ecutwfc=[],dual=[],sampling=[]):
    if ID in oneCalc['prev']:
        return oneCalc,ID
    try:
        int_samp = sorted([map(int,i.split())  for i in sampling])
        sampling=[' '.join(map(str,i)) for i in int_samp]
    except:
        pass

    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

    try:
        first_sampling = sampling[0]
    except:
        first_sampling=in_dict['K_POINTS']['__content__']
        sampling=[first_sampling]
    first_ecutwfc = ecutwfc[0]
    try:
        first_dual = dual[0]
    except:
        try:
            first_rho=int(in_dict['&system']['ecutrho'])
            first_dual=int(float(first_rho)/float(first_ecutwfc))
            dual=[first_dual]
        except:
            first_dual=4
            dual=[4]


    first_ecutrho=first_dual*first_ecutwfc
    

    #remove first entry in each list
    try:
        oneCalc['_AFLOWPI_PPTEST_ECUTWFC_']=sorted(ecutwfc)
    except:
        oneCalc['_AFLOWPI_PPTEST_ECUTWFC_']=sorted(ecutwfc)
    try:
        int_samp = sorted([map(int,i.split())  for i in sampling])
        sampling=[' '.join(map(str,i)) for i in int_samp]        
        
        oneCalc['_AFLOWPI_PPTEST_SAMPLING_']=sampling

    except:
        pass
#        int_samp = sorted([[map(int,i.split())]  for i in sampling])

    try:
        oneCalc['_AFLOWPI_PPTEST_DUAL_']=sorted(dual)
    except:
        oneCalc['_AFLOWPI_PPTEST_DUAL_']=sorted(dual)
    oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_']=[0 for x in range(len(ecutwfc))]
    oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_']=[0 for x in range(len(sampling))]
    oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_']=[0 for x in range(len(dual))]
    

    oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_'][0]=1
    oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_'][0]=1
    oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'][0]=1

    in_dict['&system']['ecutwfc']=first_ecutwfc
    in_dict['&system']['ecutrho']=first_ecutrho
    in_dict['K_POINTS']['__content__']=first_sampling

    oneCalc['_AFLOWPI_INPUT_'] = AFLOWpi.retr._joinInput(in_dict)
    oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutwfc',first_ecutwfc)
    oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutrho',first_ecutrho)
    oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'K_POINTS','__content__',first_sampling)    

    return oneCalc,ID

def _get_crawl_params(oneCalc,ID):

    new_params = AFLOWpi.retr.get_parameters(oneCalc,ID)
    param_list=['_AFLOWPI_PARAM_A_','_AFLOWPI_PARAM_B_','_AFLOWPI_PARAM_C_','_AFLOWPI_PARAM_ALPHA_','_AFLOWPI_PARAM_BETA_','_AFLOWPI_PARAM_GAMMA_',]
    param_dict={}
    for i in range(len(param_list)):
        param_dict[param_list[i]]=new_params[i]

    return param_dict


def _record__thresh(oneCalc,ID):
        crawl_params=AFLOWpi.pseudo._get_crawl_params(oneCalc,ID)
        pp_test_log_name = '_PSEUDO_TEST_LOG.log'


        in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

        ecutwfc = int(in_dict['&system']['ecutwfc'])
        ecutrho = int(in_dict['&system']['ecutrho'])
        dual = int(float(ecutrho)/float(ecutwfc))
        sampling= in_dict['K_POINTS']['__content__']
        
        ibrav=int(in_dict['&system']['ibrav'])
        
        log_params={}
        if ibrav in [1,2,3]:
            log_params={'A':crawl_params['_AFLOWPI_PARAM_A_']}


        log_params.update({'_AFLOWPI_DUAL_':dual,'_AFLOWPI_ECUTWFC_':ecutwfc,'_AFLOWPI_KPOINTS_':sampling})



        subdir=oneCalc['_AFLOWPI_FOLDER_']

        if os.path.isfile(os.path.join(subdir,'%s%s' % (ID,pp_test_log_name))):
            with open(os.path.join(subdir,'%s%s' % (ID,pp_test_log_name)),'a') as uValLog:
                uValLog.write('%s\n' % log_params)
        else:
            with open(os.path.join(subdir,'%s%s' % (ID,pp_test_log_name)),'w') as uValLog:
                uValLog.write('%s\n' % log_params)        


def _go_plot(oneCalc,ID,plot=False):
        cutoff_list=[]
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_PSEUDO_TEST_LOG.log'%ID),'r') as conv_file_obj:             
            for line in conv_file_obj.readlines():                                                         
                cutoff_list.append(ast.literal_eval(line))
        try:





            if plot==True:
                title = 'Pseudoptential Testing for %s'%AFLOWpi.retr._getStoicName(oneCalc,strip=True)
                file_name = 'PSEUDOTEST_%s_%s'%(AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID)

                AFLOWpi.pseudo.plot(cutoff_list,xaxis='_AFLOWPI_ECUTWFC_',title=title,file_name=file_name)

        except Exception,e:
            print e

        return oneCalc,ID

import time
def _check_test_conv(oneCalc,ID,thresh):

#    if oneCalc[]
#    oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_']
    #get new vals


    


    crawl_params=AFLOWpi.pseudo._get_crawl_params(oneCalc,ID)
    try:
        if '_AFLOWPI_PREV_PARAM_' not in oneCalc.keys():
            oneCalc['_AFLOWPI_PREV_PARAM_']={}



        in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
        in_dict['&control']['calculation']="'relax'"
        oneCalc['_AFLOWPI_INPUT_']=AFLOWpi.retr._joinInput(in_dict)



        num_wfc_cut = len(oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'])
        num_dual_cut = len(oneCalc['_AFLOWPI_PPTEST_DUAL_'])
        num_sampling_cut = len(oneCalc['_AFLOWPI_PPTEST_SAMPLING_'])

        if num_wfc_cut>4:
            limiter=4
        else:
            limiter=len(num_wfc_cut)
        
        p_l=['_AFLOWPI_PARAM_A_','_AFLOWPI_PARAM_B_','_AFLOWPI_PARAM_C_','_AFLOWPI_PARAM_ALPHA_','_AFLOWPI_PARAM_BETA_','_AFLOWPI_PARAM_GAMMA_',]
        conv=True


        for param in p_l:

            try:
                for i in range(limiter):


                    if numpy.abs(oneCalc['_AFLOWPI_PREV_PARAM_'][param][-(i+1)]-crawl_params[param])>thresh:
                    
                        conv=False

                        try:
                               oneCalc['_AFLOWPI_PREV_PARAM_'][param].append(crawl_params[param])
                        except:
                               oneCalc['_AFLOWPI_PREV_PARAM_'][param]=[crawl_params[param]]
                else:
                    pass
            except Exception,e:

                try:
                    oneCalc['_AFLOWPI_PREV_PARAM_'][param].append(crawl_params[param])
                except:
                    oneCalc['_AFLOWPI_PREV_PARAM_'][param]=[crawl_params[param]]
                conv=False


#       for param in p_l:
#           oneCalc['_AFLOWPI_PREV_PARAM_'][param]=crawl_params[param]
        
        oneCalc['__CRAWL_CHECK__']    
        AFLOWpi.pseudo._record__thresh(oneCalc,ID)

        printout_str='' 
        printout_str+='CONVERGED?: %s'%conv
        try:
            print AFLOWpi.run._colorize_message(printout_str,level='DEBUG',show_level=False)
        except:
            print printout_str


        if conv==True and 0 not in oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_'] and 0 not in oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_']:
            oneCalc,ID=AFLOWpi.pseudo._go_plot(oneCalc,ID,plot=True)
            in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
            try:
                pick=-1
            except:
                pick=-1
            try:
                nu_ecutwfc = oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][pick]
            except:
                nu_ecutwfc = oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][-1]
            nu_ecutwfc = oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][-1]
            
            old_ecutwfc = int(in_dict['&system']['ecutwfc'])
            old_ecutrho = int(in_dict['&system']['ecutrho'])
            nu_dual = int(float(old_ecutrho)/float(old_ecutwfc))
            nu_sampling= in_dict['K_POINTS']['__content__']
            nu_ecutrho=float(dual)*nu_ecutwfc

            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutwfc',nu_ecutwfc)
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutrho',nu_ecutrho)
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'K_POINTS','__content__',nu_sampling)
            
            printout_str='' 
            printout_str+='Converged Parameters for %s\n'%AFLOWpi.retr._getStoicName(oneCalc,strip=True)
            
            printout_str+='ECUTWFC: : %s\n' % nu_ecutwfc
            printout_str+='DUAL     : %s\n' % nu_dual
            printout_str+='SAMPLING : %s\n' % nu_sampling
            try:
                print AFLOWpi.run._colorize_message(printout_str,level='DEBUG',show_level=False)
            except:
                print printout_str
            print 
            

            return oneCalc,ID

        elif 0 not in oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_'] and 0 not in oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_'] and 0 not in oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'] and conv == False:
            logging.info('Pseudotesting convergence not achieved. Exiting')
            raise SystemExit
        else:
            #if converged but dual and sampling to check still shift those and keep ecutwfc the same
            
            if conv==True:
                
                if 0 in oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_']:
                    shift_type='dual'
                
                elif 0 in oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_']:
                    shift_type='sampling'
                else:
                    print 'this should not happen'
            else:
                shift_type='ecutwfc'

            oneCalc,ID=AFLOWpi.pseudo._shift_cutoffs(oneCalc,ID,shift_type=shift_type)

            try:
                oneCalc['__execCounter__']=0
                try:
                    del oneCalc['__CRAWL_CHECK__']
                    try:
                        del oneCalc['__gridMin_iteration__']
                    except:
                        pass
                except Exception,e:
                    pass
            except Exception,e:
                pass
            return oneCalc,ID
    except Exception,e:

        return oneCalc,ID

def _shift_cutoffs(oneCalc,ID,shift_type='wfc'):
    print 

    oneCalc['_AFLOWPI_PREV_SHIFT_TYPE']=shift_type
    

    printout_str=''
    if shift_type=='ecutwfc':
        try:

            for i in range(len(oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'])):

                if oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'][i]==0:
                    
                    oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'][i]=1
                    new_ecutwfc=oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][i]
                    break
                else:
                    new_ecutwfc=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ecutwfc']

            ecutw = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ecutwfc']
            ecutr = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ecutrho']
            
            dual= float(ecutr)/float(ecutw)

            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutwfc',new_ecutwfc)    
            new_ecutr=int(new_ecutwfc*dual)

            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutrho',new_ecutr)    
            printout_str+='Shifting wavefunction energy cutoff\n'
            
        except Exception,e:
            print e
            raise SystemExit
    if shift_type=='dual':
        try:
            
            oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'] = [0 for i in oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_']]
            oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'][0]=1

            
            for i in range(len(oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_'])):

                if oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_'][i]==0:
                    oneCalc['_AFLOWPI_PPTEST_DUAL_TRACKER_'][i]=1
                    new_dual=oneCalc['_AFLOWPI_PPTEST_DUAL_'][i]
                    INDEX=i
                    break

#            oneCalc['_AFLOWPI_PPTEST_ECUTWFC_']=[x for x in reversed(oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][:INDEX+1])]

            oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'] = [0 for i in oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_']]
            oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'][0]=1
            del oneCalc['_AFLOWPI_PREV_PARAM_']

            new_ecutwfc=oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][0]
            
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutwfc',new_ecutwfc)    
#            new_ecutr=int(new_ecutwfc*dual)
            new_ecutrho=int(new_ecutwfc*float(new_dual))


            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutrho',new_ecutrho)

            printout_str+='Shifting dual\n'


            
        except Exception,e:

            print e
            raise SystemExit
            pass 
    if shift_type=='sampling':

        try:
            for i in range(len(oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_'])):

                if oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_'][i]==0:
                    oneCalc['_AFLOWPI_PPTEST_SAMPLING_TRACKER_'][i]=1
                    new_sampling=oneCalc['_AFLOWPI_PPTEST_SAMPLING_'][i]
                    INDEX=i
                    break
                #else:
                #    return oneCalc,ID
                #    new_sampling=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['K_POINTS']['__content__']
        except:

            print e
            raise SystemExit
        


#        oneCalc['_AFLOWPI_PPTEST_ECUTWFC_']=[x for x in reversed(oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][:INDEX+1])]
        
        oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'] = [0 for i in oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_']]
        oneCalc['_AFLOWPI_PPTEST_ECUTWFC_TRACKER_'][0]=1
        
        del oneCalc['_AFLOWPI_PREV_PARAM_']
        
        ecutwfc=oneCalc['_AFLOWPI_PPTEST_ECUTWFC_'][0]
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'K_POINTS','__content__',new_sampling)    
        ecutw = float(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ecutwfc'])
        ecutr = float(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ecutrho'])
        dual=int(ecutr/ecutw)
        printout_str+='Shifting sampling\n'

        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutwfc',ecutwfc)    
        new_ecutrho=int(ecutwfc*dual)
        oneCalc['_AFLOWPI_PREV_PARAM_']
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','ecutrho',new_ecutrho)    


    try:
        in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

        nu_ecutwfc = int(in_dict['&system']['ecutwfc'])
        nu_ecutrho = int(in_dict['&system']['ecutrho'])
        nu_dual = int(float(nu_ecutrho)/float(nu_ecutwfc))
        nu_sampling= in_dict['K_POINTS']['__content__']


        printout_str+='New testing Parameters for %s\n'%AFLOWpi.retr._getStoicName(oneCalc,strip=True)
        printout_str+='ECUTWFC: : %s\n' % nu_ecutwfc
        printout_str+='DUAL     : %s\n' % nu_dual
        printout_str+='SAMPLING : %s\n' % nu_sampling
        try:
            print AFLOWpi.run._colorize_message(printout_str,level='DEBUG',show_level=False)
        except:
            print printout_str


    except Exception,e:
        print e
    
    return oneCalc,ID



def _getCutOffs(calcs):
    '''
    Looks through a random 'sample' calculations in the dictionary of dictionaries of calculations
    to see which cutoff variables are defined in the set of calculations and returns a list of
    those cutoffs

    Arguments:
     - calcs -- Dictionary of dictionaries of calculations

    '''
    sampleDict = calcs[random.choice(calcs.keys())]
    key=[]
    try:
	    if '_AFLOWPI_KPOINTS_' in sampleDict:
		key.append('_AFLOWPI_KPOINTS_')
	    if '_AFLOWPI_ECUTR_'  in sampleDict: 
		if '_AFLOWPI_DUAL_' not in sampleDict:
			key.append('_AFLOWPI_ECUTR_')
		else:
			key.append('_AFLOWPI_DUAL_')
	    if '_AFLOWPI_ECUTW_' in sampleDict:
		key.append('_AFLOWPI_ECUTW_')

	    if len(key)==0:
		    return []
	    return key
    except:
	    return []

def _splitCalcs(calcs,splitVars=''):
	'''
	splits a set of calculations by keyword and its respective values and
	returns the split set of calculations as a list of the split set.

	Arguments:
	      calcs (dict): dictionary of dictionary of calculations that is to be split

	Keyword Arguments:
	      splitVars (list): a list or tuple of strings of AFLOWpi variable(s) in the
                                set of calculations that you want to split the set on

        Returns:
              splitCalcList (list): a list of dictionaries of dictionaries

	'''

	if splitVars=='' or len(splitVars)==0 or splitVars==['']:
		return [calcs]

	key=[]
	value=[]
	calcsCopy = copy.deepcopy(calcs)
	sampleDict = calcsCopy[random.choice(calcsCopy.keys())]   
	if sampleDict.has_key('_AFLOWPI_DUAL_'):
		try:
			splitVars.remove('_AFLOWPI_ECUTR_')
			splitVars.append('_AFLOWPI_DUAL_')
		except ValueError:
			pass

	for var in splitVars:
	    if var in sampleDict:
		inAllCalcs = []
		for ID,oneCalc in calcsCopy.iteritems(): #cycles through all the calculations
		    inAllCalcs.append(oneCalc[var]) #makes a list of all the values for all your calcs
		key.append(var)

		listCalcs = list(set(inAllCalcs))

		sortedListCalcs= sorted(listCalcs)
		value.append(tuple(sortedListCalcs)) #a list tuples containing the unique values of the variables you are splitting on in your calcs 


	keyTuple   = tuple(key)
	valueTuple = tuple(value)
	cutoffTuple =  tuple(it.product(*valueTuple))
	constrList=[]
	for i in cutoffTuple:
	    constraintDict=OrderedDict()
	    for k,v in enumerate(i):
		key = keyTuple[k]
		constraintDict[key] = v
	    constrList.append(constraintDict)

	splitCalcList=[]
	for constraint in constrList:
	    splitOn =OrderedDict()

	    for ID,oneCalc in calcsCopy.iteritems():
		if len([i for i in constraint if constraint[i]==oneCalc[i]]) == len(constraint): #fancy way of looping 'if' statements
		    splitOn.update({ID:oneCalc})

	    splitOnCopy=copy.deepcopy(splitOn)
	    splitCalcList.append(splitOnCopy)

	return splitCalcList

	    

def _grabEnergyOut(calcs):
    """
    <CONSIDER MOVING TO retr.py>
    Goes in every subdirectory of the calculation and searches for the final energy of the calculation
    and returns a new copy of the input dictionary that includes the final energy.

    Arguments:
          calcs (dict): Dictionary of Dictionaries of the calculations.

    Keyword Arguments:
          None

    MASTER_ENERGY (dict): Dictionary for the calc set with the energy added to each calc in the set
    """
    
    MASTER_ENERGY = copy.deepcopy(calcs)
    energyRegex = re.compile(r'(?:(?:(?:(?:\!\s+)total)|(?:Final)) en\w+\s*=\s+(.+?)Ry)',re.MULTILINE)
    
    for ID,oneCalc in calcs.iteritems():
        try:

		if os.path.exists(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)):
		    with file(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID),'r') as outFile:
			outFileString=outFile.read()
			#searches for either total energy or final energy in the scf output and adds it to the output dictionary
			finalEnergy=energyRegex.findall(outFileString)

			if len(finalEnergy):
			    energyFloat  = float(finalEnergy[-1])

			    MASTER_ENERGY[ID]['Energy']= energyFloat
			else: #if the energy can not be found the test entry is deleted from the output dictionary
			    outCalcPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)
			    logging.warning('could not get energy. check output file: %s' % outCalcPath)
			    print 'could not get energy. check output file: %s' % outCalcPath
                            MASTER_ENERGY[ID]['Energy']=None
                            calcs[ID]['Energy']=0.0
			    del MASTER_ENERGY[ID] 
#			    del calcs[ID]
		else:
                    MASTER_ENERGY[ID]['Energy']=None
		    del MASTER_ENERGY[ID]            
	except:
		outCalcPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)
		logging.warning('could not get energy. check output file: %s' % outCalcPath)
		print 'could not get energy. check output file: %s' % outCalcPath
    return MASTER_ENERGY





#########################################################################################################
try:
	import scipy
except Exception,e:
	try:
		from scipy.optimize import differential_evolution
	except Exception,e:
		try:
			from scipy.optimize import *
		except Exception,e:
			pass

from scipy.optimize import minimize
from numpy import ones
from numpy import array
from numpy import zeros
from numpy import triu_indices
from numpy import diag_indices
from numpy import dot
#########################################################################################################

from scipy.interpolate import griddata
from scipy.ndimage import map_coordinates,spline_filter
def _cubicInterpolation(variable_array,solution_array):

	'''try to convert a numpy array into a list if
	this fails assume that the inputs are lists'''
	try:
		variable_array=variable_array.tolist()
		solution_array=solution_array.tolist()
	except Exception,e:
		logging.info('could not convert variable_array or solution_array into a list in AFLOWpi.pseudo._cubicInterpolation. Possibly they are already lists?')

	'''make a unique sorted list for our grid from the data points'''
	sortedVariableArray=[sorted(list(set(x))) for x in variable_array]
	insidePad=[]
	'''
	for every entry in the range of a variable 
	(excluding the end and beginning of the range)
	insert a value in between them that's the average of the two.)
	'''
	for var in range(len(sortedVariableArray)):
		tempPad=[]
		for x in range(len(sortedVariableArray[var])):
			even = sortedVariableArray[var][x]			
			try:
				odd=(sortedVariableArray[var][x]+sortedVariableArray[var][x+1])/2.0
				tempPad.append(even)
				tempPad.append(odd)
			except:
				tempPad.append(even)
		
		insidePad.append(sorted(tempPad))
	'''cartesian product to recreate our new M*(M[N]*2-1) long list'''
	variableListPadded = list(it.product(*insidePad))
	paddedGrid = list(it.product(*insidePad))
	mapCoordVars= zip(*paddedGrid)
 	paddedGrid =  numpy.asarray(paddedGrid)
	variable_array = numpy.asarray(variable_array).T
	'''mask bad energy values from the calculations'''
	energyMatrix = numpy.asarray(solution_array)	
	energyMatrix = numpy.ma.masked_values(energyMatrix,0.0)

#	reshapeDataParam=numpy.array([len(x) for x in sortedVariableArray])
#	reshapeDataParam2=numpy.array([len(x)*2-1 for x in sortedVariableArray])
#	energyMatrixReshape=numpy.reshape(energyMatrix,reshapeDataParam)

	
	grid_z0 = griddata(variable_array, energyMatrix,paddedGrid, method='linear')
        
	return paddedGrid.T.tolist(),grid_z0.tolist()


import scipy
def _cubicInterpolation_onePoint(point,variable_array,solution_array):

	'''try to convert a numpy array into a list if
	this fails assume that the inputs are lists'''
	try:
		variable_array=variable_array.tolist()
		solution_array=solution_array.tolist()
	except Exception,e:
		logging.info('could not convert variable_array or solution_array into a list in AFLOWpi.pseudo._cubicInterpolation. Possibly they are already lists?')

	'''make a unique sorted list for our grid from the data points'''
	sortedVariableArray=[sorted(list(set(x))) for x in variable_array]
        insidePad=[]
	'''
	for every entry in the range of a variable 
	(excluding the end and beginning of the range)
	insert a value in between them that's the average of the two.)
	'''

	'''cartesian product to recreate our new M*(M[N]*2-1) long list'''
	try:
		energyMatrix = numpy.asarray(solution_array)	
		energyMatrix = numpy.ma.masked_values(energyMatrix,0.0)



		'''round all the values of the variable_array to account for precision errors'''

		variable_array = numpy.asarray(variable_array)

		sortedVariableArray=[sorted(list(set(x))) for x in variable_array]
#		reshapeDataParam=#
#		energyMatrixReshape=numpy.reshape(energyMatrix,reshapeDataParam).T



#		variable_array=numpy.around(variable_array,decimals=6)


#		print reshapeDataParam


#		energyMatrix=numpy.reshape(energyMatrixReshape,len(energyMatrix))
	#	gg=numpy.array((1,),dtype=numpy.float64)
	#	gg=numpy.zeros(1, dtype=float, order='C')






		gg=[]
	# 	for x in range(len(variable_array[0])):
	# 		gg.append([variable_array[0][x],energyMatrixReshape[x]])
	# 	energyMatrixReshape=numpy.asarray(gg)

	# 	print energyMatrixReshape.T
	#	print point

		point = numpy.array(point)


#		print energyMatrix


		if variable_array.shape[1]!=1:
			point.T
#			energyMatrix=numpy.reshape(energyMatrix,(1,len(energyMatrix)))

		if variable_array.shape[1]==1:
			variable_array=variable_array.T[0]
#		print variable_array
#		print  variable_array.shape

#		rbfi = scipy.interpolate.Rbf(variable_array,energyMatrix,function='multiquadric')
	#	print rbfi
	#	print 
	#	grid_z0 = map_coordinates(energyMatrixReshape, [point],order=5)
	#	print grid_z0
	#	grid_z0 = (variable_array, energyMatrix,point, method='linear')


		if len(variable_array.shape)<3:
			method='cubic'
                        
		else:
			method='linear'


                if variable_array.shape<2:
                    f2 = interp1d(variable_array,solution_array, kind='linear')
                    sol = f2(point)

                elif variable_array.shape==2:
                    f2 = interp2d(variable_array,solution_array, kind='cubic')
                    sol = f2(point)
                else:
		    sol = griddata(variable_array, energyMatrix,point, method=method)


#		sol = rbi(point)

		try:
                        
			return sol[0][0]
		except:
			try:
				return sol[0]
			except:
				return sol

	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
		raise SystemExit

def _getMatrices(calcs,fitVars=None,options=None,Energy=True):
#    print 'entering __getMatrices'
    logging.debug('entering __getMatrices')

    totalList = []
    energyList = [] # pulls energy values from the dictionary of dictionaries of calculations and creates a list of them
    orderedMaster = OrderedDict(calcs)
    if Energy==True:
        for ID,oneCalc in orderedMaster.iteritems():
            energyList.append(oneCalc['Energy'])
    for variables in fitVars:
        varValueList = []
        for ID,oneCalc in orderedMaster.iteritems():
            varValueList.append(float(oneCalc[variables]))
	totalList.append(varValueList)
        
    
    
    energyMatrix = array(energyList)
    

    variableList = array(zip(*totalList))

#    print 'exiting __getMatrices'
    logging.debug('exiting __getMatrices')
    if Energy==True:
        return energyMatrix,variableList
    else:
        return variableList

def _getMin(calcs,fitVars=None,options=None,bulk_modulus=False,minimize_var="Energy"):
    '''
    calculates the minimum energy value from the fit variables you supply. 
    Only picks out the energy values that satisfy the cutoff/kpoint values 
    for one combination of them that you are testing 
    (i.e. would calculate the minimum energy from calculations that satisfy 
    {_AFLOWPI_KPOINTS_ = '4 4 4 1 1 1', _AFLOWPI_ECUTW_ = 20}  

    Arguments:
          calcs (dict): Dictionary of Dictionaries of the calculations.
    Keyword Arguments
          fitVars (list): a tuple of the independent variables you want to make the 
                          fit and find the min energy for
          options (dict): additional options passed to L-BFGS-B minimization 
                          (see scipy documentation for scipy.optimize.minimize 
                          for more details)

    Returns:
         minEnergy (float): the minimum energy of the fit
         outDict (dict): dictionary containing fit information

    '''

    def fitFunction(x,*p):     
	    if not x.size:

		    return 
	    try:
		n = x.size



		f0 = p[0]
		x0 = p[1:n+1]
		M = zeros((n,n),dtype=numpy.float64)
		D = zeros((n,n),dtype=numpy.float64)
		iu = triu_indices(n)
		M[iu] = p[n+1:]
		di = diag_indices(n)
		D[di] = M[di]

		M = M-D + M.T
		l = lambda z:  f0+dot(z-x0,dot(M,z-x0))
		# print di
		# print D
		# print M
		# print array(l(x))
		# print f0
		# print p
		# print n
		# print x0
		# print M
		# print D
		# print iu
		# print M
		return array(l(x))
	    except Exception,e:
		    AFLOWpi.run._fancy_error_log(e)
		    return 


    '''
    generic form of the taylor expansion for n-dim 
    '''
    def fitForm(x,*p):     
        n = x[0].size
        f0 = numpy.array(p[0],dtype=numpy.float64)
        x0 = array(p[1:n+1])

        M = zeros((n,n),dtype=numpy.float64)
        D = zeros((n,n),dtype=numpy.float64)
        iu = triu_indices(n)
        M[iu] = p[n+1:]
        di = diag_indices(n)
        D[di] = M[di]
        M = M-D + M.T
        l = lambda z:  f0 +dot(z-x0,dot(M,z-x0))
        return array(map(l,x))


    logging.debug('entering __getMin')

    totalList = []
    energyList = [] # pulls energy values from the dictionary of dictionaries of calculations and creates a list of them
    for ID,oneCalc in calcs.iteritems():
        try:
            oneCalc['Energy']
        except:
            calcs = AFLOWpi.pseudo._grabEnergyOut(calcs)
            break
    orderedMaster = OrderedDict(calcs)
    for ID,oneCalc in orderedMaster.iteritems():
        energyList.append(oneCalc[minimize_var])
    

#    for variables in fitVars: #cycling through the independent variables you are fitting with to prep the data for the fit
        
     
    '''
    cycles through all the calculations in your dictionary of dictionaries
    and checks if the calculation satisfies the cutoff/kpoint combination
    that you want to find the minimum for. It preps the data for the 
    minimum energy calculation
    '''


    for variables in fitVars:
        varValueList = []
        for ID,oneCalc in orderedMaster.iteritems():
            varValueList.append(float(oneCalc[variables]))
	totalList.append(varValueList)
    
    try: 
	    
#	    totalList,energyList = __cubicInterpolation_onePoint(totalList,energyList,x)
	    variableList = numpy.asarray(totalList).T


	    energyMatrix = numpy.asarray(energyList)


	    '''mask the bad energy values from failed calcs'''
	    energyMatrix = numpy.ma.masked_values(energyMatrix,0.0)
	    #chooses the median value for each of the fitting variables as starting points for the minimization



    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)


    

    '''
     fits a curve which is a n-dim second order taylor expansion to the gridded data set
    '''
    n = variableList[0].size
    m = 1 + n+ n*(n+1)/2
    guess = ones(m) 

    '''
    gets bounds of the data set for the bounded BFGS minimization and then
    finds the minimum value of the fitted function
    '''
    ystart = numpy.array([numpy.float64(numpy.median(variable)) for variable in variableList.T])
    
    fitBounds = list(tuple([numpy.amin(variable)+0.00000001,numpy.amax(variable)-0.00000001]) for variable in variableList.T)
    try:
	    pop, cov = curve_fit(fitForm,variableList,energyMatrix,guess,maxfev=10000)
    except:
            
	    try:
		    pop, cov = curve_fit(fitForm,variableList,energyMatrix,guess,maxfev=100000)
	    except:
		    try:
			    pop, cov = curve_fit(fitForm,variableList,energyMatrix,guess,maxfev=1000000)
		    except Exception,e:
			    AFLOWpi.run._fancy_error_log(e)


    try:
	    pop=tuple(pop.tolist())
    except Exception,e: 
	    AFLOWpi.run._fancy_error_log(e)



    try:
        logging.info('minimizing using differential evolution')

        fitMin = scipy.optimize.differential_evolution(AFLOWpi.pseudo._cubicInterpolation_onePoint,fitBounds,popsize=1000,
                                                       args=(variableList,energyMatrix),disp=False,tol=0.000001,
                                                       maxiter=100)
        
#

  #      fitMin = minimize(fitFunction,ystart,args=(pop),method='L-BFGS-B',bounds=fitBounds,tol=0.00001,options={'disp':False})
    except Exception,e: 
#            AFLOWpi.run._fancy_error_log(e)
            logging.info('minimizing using bfgs')
	    AFLOWpi.run._fancy_error_log(e)

	    try:
                fitMin = scipy.optimize.differential_evolution(fitFunction,fitBounds,popsize=500,args=(pop),
                                                               disp=True,tol=0.00000001,maxiter=100) 	    
                                                               

	    except Exception,e: 
                    AFLOWpi.run._fancy_error_log(e)

		    try:
                        fitMin = scipy.optimize.minimize(AFLOWpi.pseudo._cubicInterpolation_onePoint,ystart,
                                                         args=(variableList,energyMatrix),method='L-BFGS-B',
                                                         bounds=fitBounds,options={'disp':False,
                                                         'gtol':0.0001}) 
		    except Exception,e:
			    AFLOWpi.run._fancy_error_log(e)

			    try:
                                fitMin = scipy.optimize.minimize(AFLOWpi.pseudo._cubicInterpolation_onePoint,ystart,
                                                                 args=(variableList,energyMatrix),method='L-BFGS-B',
                                                                 bounds=fitBounds,options={'disp':False,
                                                                                           'gtol':0.00001,
                                                                                           'eps':0.000000001},) 
			    except Exception,e:
                                AFLOWpi.run._fancy_error_log(e)
                                pass


    if bulk_modulus==True:
	    try:
		    with open('../AFLOWpi/bulk_modulus.dat','w') as bulk_mod_file:
			    bulk_mod_file.write('%s'%repr(pop))
			    bulk_mod_file.write('\n')
	    except Exception,e:
		    AFLOWpi.run._fancy_error_log(e)
		    

#    print fitMin1
#    if options!=None:
#	        fitMin = minimize(fitFunction,ystart,args=(pop),method='L-BFGS-B',bounds=fitBounds,options=options,disp=True)
 #   else:


    minEnergy=fitMin.fun

    blankDict = {}
    outDict=OrderedDict(blankDict)
    '''
    puts the results from the fitting and minimization in a dictionary
    and returns it to be further processed
    '''

    logging.debug('exiting __getMin')

    for item in range(len(fitMin.x)):
        outDict[fitVars[item]]=fitMin.x[item]




    return minEnergy,outDict



def  getMinimization(origCalcs,fitVars=None,options=None,runlocal=False,faultTolerant=True,minimize_var="Energy"):
	if runlocal==True:
		return AFLOWpi.pseudo._getMinimization(origCalcs,fitVars=fitVars,options=options,minimize_var=minimize_var)
	##########################################################################################
	###NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED ###
	##########################################################################################
	else:
		command = '''
         center = AFLOWpi.pseudo._getCenter(calcs,fitVars=%s,options=%s)
         

         AFLOWpi.pseudo.plot(resultList,xaxis='',xtitle=None,ytitle=None,title=None,rename=None)

''' % (fitVars,options)

		AFLOWpi.prep.runAfterAllDone(origCalcs,command,faultTolerant=faultTolerant)
	
	##########################################################################################
	###NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED NOT COMPLETED ###
	##########################################################################################


def  crawlingMinimization(calcs,options=None,faultTolerant=True,constraint=None,thresh=0.001,initial_variance=0.05,grid_density=10,mult_jobs=False,final_minimization='relax',calc_type='relax'):


        constraintString='None'
        if constraint==None:
            constraintString='None'
        else:
            try:
                if type(constraint[0][1])==type('string'):


                
                    constraintString='('
                    constraintString+='("%s","%s"),' % (constraint[0],constraint[1])
                    constraintString+=')'


            except Exception,e:

                    print 'constraint list not formatted properly. must be in form constraint=(("constraint1","parameter_constrained1"),("constraint2","parameter_constrained2"),) ....exiting'
                    logging.warning('constraint list not formatted properly. must be in form constraint=(("constraint1","parameter_constrained1"),("constraint2","parameter_constrained2"),) ....exiting')
                    AFLOWpi.run._fancy_error_log(e)
                    raise SystemExit

	fitVars=None
#########$$$$$$$$$$$$$$$$$$
        if calc_type.lower() != 'relax' and calc_type.lower() != 'scf':
            print 'only "scf" and "relax" supported for calc_type'
            sys.exit(0)
        AFLOWpi.prep._addToAll(calcs,'PREPROCESSING',"""oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation',"%s")""" %repr(calc_type))


        add_outfile_task="""AFLOWpi.prep._addToAll(calc_subset,'PREPROCESSING',"outFile=os.path.join('%s','_%s.py')"%(oneCalc['_AFLOWPI_FOLDER_'],ID))"""
        add_vars='AFLOWpi.pseudo._crawl_min_vars_to_calc(calc_subset)'
        run_sub='AFLOWpi.run.scf(calc_subset,exit_on_error=False)'

        task_list=[add_vars,add_outfile_task,run_sub]

        for ID,oneCalc in calcs.iteritems():
            refFileDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
            paramList=[]
            if '_AFLOWPI_ORIG_CDM1_' in oneCalc.keys():
                    paramList.append('A')
            if '_AFLOWPI_ORIG_CDM2_' in oneCalc.keys():
                    paramList.append('B')
            if '_AFLOWPI_ORIG_CDM3_' in oneCalc.keys():
                    paramList.append('C')
            if '_AFLOWPI_ORIG_CDM4_' in oneCalc.keys():
                    paramList.append('alpha')
            if '_AFLOWPI_ORIG_CDM5_' in oneCalc.keys():
                    paramList.append('beta')
            if '_AFLOWPI_ORIG_CDM6_' in oneCalc.keys():
                    paramList.append('gamma')
            break

        fitVars=['_AFLOWPI_PARAM_%s_'%param for param in paramList]

        shift_check="""AFLOWpi.pseudo._shiftGrid(calcs,outFile,fitVars=%s,options=%s,constraint=%s,thresh=%s,mult_jobs=%s)"""%(fitVars,options,constraint,thresh,mult_jobs)
        gen_set_command = 'AFLOWpi.pseudo._crawlingMinimization(oneCalc,ID,faultTolerant=%s,constraint=%s,thresh=%s,initial_variance=%s,steps=%s,mult_jobs=%s,)'%(faultTolerant,constraint,thresh,initial_variance,grid_density,mult_jobs,)

        calcs=AFLOWpi.prep.prep_split_step(calcs,gen_set_command,subset_tasks=task_list,mult_jobs=mult_jobs,substep_name='GRID_MIN',keep_file_names=False,clean_input=False,check_function=shift_check,fault_tolerant=True)


	if final_minimization not in ['vc-relax','relax'] and final_minimization != None:
		final_minimization='scf'


        if final_minimization!=None:
#            AFLOWpi.run._skeletonRun(calcs)
            AFLOWpi.prep._addToAll(calcs,'RUN',"oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','%s')" %final_minimization)
            for ID,oneCalc in calcs.iteritems():
                oneCalc['__execCounterBkgrd__']+=1
            AFLOWpi.run.scf(calcs)

	return calcs




###############################################################################################################

###############################################################################################################

########################################################################################################################################################################################################################################


def _crawlingMinimization(oneCalc,ID,fitVars=None,options=None,faultTolerant=True,initial_variance=0.15,steps=10,constraint=None,thresh=0.001,mult_jobs=True):


	chain_index=1

	chain_logname='step_%02d'%1
        checkBool=False
	if not checkBool:
		assocInputVars=OrderedDict()
		try:
			refFileDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
		except:
			pass

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


		centValList=[]
		paramList=[]
                if '_AFLOWPI_ORIG_CDM1_' in oneCalc.keys():
                    paramList.append('A')
                    centValList.append(float(oneCalc['_AFLOWPI_ORIG_CDM1_']))
                if '_AFLOWPI_ORIG_CDM2_' in oneCalc.keys():
                    paramList.append('B')
                    centValList.append(float(oneCalc['_AFLOWPI_ORIG_CDM2_']))
                if '_AFLOWPI_ORIG_CDM3_' in oneCalc.keys():
                    paramList.append('C')
                    centValList.append(float(oneCalc['_AFLOWPI_ORIG_CDM3_']))
                if '_AFLOWPI_ORIG_CDM4_' in oneCalc.keys():
                    paramList.append('alpha')
                    centValList.append(float(oneCalc['_AFLOWPI_ORIG_CDM4_']))
                if '_AFLOWPI_ORIG_CDM5_' in oneCalc.keys():
                    paramList.append('beta')
                    centValList.append(float(oneCalc['_AFLOWPI_ORIG_CDM5_']))
                if '_AFLOWPI_ORIG_CDM6_' in oneCalc.keys():
                    paramList.append('gamma')
                    centValList.append(float(oneCalc['_AFLOWPI_ORIG_CDM6_']))

		# if 'celldm(1)' in refFileDict['&system'].keys():
		# 	paramList.append('A')
		# 	centValList.append(refFileDict['&system']['celldm(1)'])
		# if 'celldm(2)' in refFileDict['&system'].keys():
		# 	paramList.append('B')
		# 	centValList.append(refFileDict['&system']['celldm(2)'])
		# if 'celldm(3)' in refFileDict['&system'].keys():
		# 	paramList.append('C')
		# 	centValList.append(refFileDict['&system']['celldm(3)'])
		# if 'celldm(4)' in refFileDict['&system'].keys():
		# 	paramList.append('alpha')
		# 	centValList.append(refFileDict['&system']['celldm(4)'])
		# if 'celldm(5)' in refFileDict['&system'].keys():
		# 	paramList.append('beta')
		# 	centValList.append(refFileDict['&system']['celldm(5)'])
		# if 'celldm(6)' in refFileDict['&system'].keys():
		# 	paramList.append('gamma')
		# 	centValList.append(refFileDict['&system']['celldm(6)'])



		stepList=[steps for x in range(len(paramList))]
		amountList=[initial_variance for x in range(len(paramList))]



		variedCalcs=AFLOWpi.prep.varyCellParams(oneCalc,ID,param=paramList,amount=amountList,steps=stepList,constraint=constraint)
                return variedCalcs
###########################################################################################################################
def _crawl_min_vars_to_calc(calc_subset):
    for ID,oneCalc in calc_subset.iteritems():
        refFileDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
        centValList=[]
        paramList=[]
        if 'celldm(1)' in refFileDict['&system'].keys():
                paramList.append('A')
                centValList.append(refFileDict['&system']['celldm(1)'])
        if 'celldm(2)' in refFileDict['&system'].keys():
                paramList.append('B')
                centValList.append(refFileDict['&system']['celldm(2)'])
        if 'celldm(3)' in refFileDict['&system'].keys():
                paramList.append('C')
                centValList.append(refFileDict['&system']['celldm(3)'])
        if 'celldm(4)' in refFileDict['&system'].keys():
                paramList.append('alpha')
                centValList.append(refFileDict['&system']['celldm(4)'])
        if 'celldm(5)' in refFileDict['&system'].keys():
                paramList.append('beta')
                centValList.append(refFileDict['&system']['celldm(5)'])
        if 'celldm(6)' in refFileDict['&system'].keys():
                paramList.append('gamma')
                centValList.append(refFileDict['&system']['celldm(6)'])




        cdmDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']
        refDictDict = AFLOWpi.retr._splitInput(oneCalc['__refFile__'])
        for param in paramList:

            if param=='A':
                val=float(cdmDict['celldm(1)'])
                oneCalc['_AFLOWPI_PARAM_%s_'%param]=val
                refDictDict['&system']['celldm(1)']='_AFLOWPI_PARAM_%s_'%param
            try:
                if param=='B':
                    val = float(cdmDict['celldm(1)'])*float(cdmDict['celldm(2)'])
                    oneCalc['_AFLOWPI_PARAM_%s_'%param]=val
                    refDictDict['&system']['celldm(2)']='_AFLOWPI_PARAM_%s_'%param
            except:
                pass
            try:
                if param=='C':
                    val = float(cdmDict['celldm(1)'])*float(cdmDict['celldm(3)'])
                    oneCalc['_AFLOWPI_PARAM_%s_'%param]=val
                    refDictDict['&system']['celldm(3)']='_AFLOWPI_PARAM_%s_'%param
            except:
                pass
            try:
                if param=='alpha':
                    val = numpy.arccos(float(cdmDict['celldm(4)']))*(180/numpy.pi)
                    oneCalc['_AFLOWPI_PARAM_%s_'%param]=val
                    refDictDict['&system']['celldm(4)']='_AFLOWPI_PARAM_%s_'%param
            except:
                pass
            try:
                if param=='beta':
                    val = numpy.arccos(float(cdmDict['celldm(5)']))*(180/numpy.pi)
                    oneCalc['_AFLOWPI_PARAM_%s_'%param]=val
                    refDictDict['&system']['celldm(5)']='_AFLOWPI_PARAM_%s_'%param
            except:
                pass
            try:
                if param=='gamma':
                    val = numpy.arccos(float(cdmDict['celldm(6)']))*(180/numpy.pi)
                    oneCalc['_AFLOWPI_PARAM_%s_'%param]=val
                    refDictDict['&system']['celldm(6)']='_AFLOWPI_PARAM_%s_'%param
            except:
                pass



        oneCalc['__refFile__']=AFLOWpi.retr._joinInput(refDictDict)
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)


    return calc_subset
###########################################################################################################################


def _add_shifter(calc_subset,fitVars=None,options=None,faultTolerant=True,initial_variance=0.15,steps=10,constraint=None,thresh=0.001,mult_jobs=True):
    for ID_new,oneCalc_new in calc_subset.iteritems(): 
    
                cdmDict = AFLOWpi.retr._splitInput(oneCalc_new['_AFLOWPI_INPUT_'])['&system']
                refDictDict = AFLOWpi.retr._splitInput(oneCalc_new['__refFile__'])


                paramList=[]
                if 'celldm(1)' in refFileDict['&system'].keys():
                        paramList.append('A')
                if 'celldm(2)' in refFileDict['&system'].keys():
                        paramList.append('B')
                if 'celldm(3)' in refFileDict['&system'].keys():
                        paramList.append('C')
                if 'celldm(4)' in refFileDict['&system'].keys():
                        paramList.append('alpha')
                if 'celldm(5)' in refFileDict['&system'].keys():
                        paramList.append('beta')
                if 'celldm(6)' in refFileDict['&system'].keys():
                        paramList.append('gamma')


                fitVars=[]
                for param in paramList:
                    fitVars.append('_AFLOWPI_PARAM_%s_'%param)


                for param in paramList:
                        if param=='A':
                                val=float(cdmDict['celldm(1)'])
                                oneCalc_new['_AFLOWPI_PARAM_%s_'%param]=val
                                refDictDict['&system']['celldm(1)']='_AFLOWPI_PARAM_%s_'%param
                        if param=='B':
                                val = float(cdmDict['celldm(1)'])*float(cdmDict['celldm(2)'])
                                oneCalc_new['_AFLOWPI_PARAM_%s_'%param]=val
                                refDictDict['&system']['celldm(2)']='_AFLOWPI_PARAM_%s_'%param
                        if param=='C':
                                val = float(cdmDict['celldm(1)'])*float(cdmDict['celldm(3)'])
                                oneCalc_new['_AFLOWPI_PARAM_%s_'%param]=val
                                refDictDict['&system']['celldm(3)']='_AFLOWPI_PARAM_%s_'%param
                        if param=='alpha':
                                val = numpy.arccos(float(cdmDict['celldm(4)']))*(180/numpy.pi)
                                oneCalc_new['_AFLOWPI_PARAM_%s_'%param]=val
                                refDictDict['&system']['celldm(4)']='_AFLOWPI_PARAM_%s_'%param
                        if param=='beta':
                                val = numpy.arccos(float(cdmDict['celldm(5)']))*(180/numpy.pi)
                                oneCalc_new['_AFLOWPI_PARAM_%s_'%param]=val
                                refDictDict['&system']['celldm(5)']='_AFLOWPI_PARAM_%s_'%param
                        if param=='gamma':
                                val = numpy.arccos(float(cdmDict['celldm(6)']))*(180/numpy.pi)
                                oneCalc_new['_AFLOWPI_PARAM_%s_'%param]=val
                                refDictDict['&system']['celldm(6)']='_AFLOWPI_PARAM_%s_'%param






                oneCalc_new['__refFile__']=AFLOWpi.retr._joinInput(refDictDict)
                AFLOWpi.prep._saveOneCalc(oneCalc_new,ID_new)

		AFLOWpi.prep.updatelogs(variedCalcs,logname=chain_logname)
		try:
			if constraint==None:
				constraintString='None'
			else:
				constraintString='('
				for constraints in range(len(constraint_type_list)):
					constraintString+='("%s","%s"),' % (constraint_type_list[constraints],constraint_var_list[constraints])
				constraintString+=')'

		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)
			constraintString='None'




#		We don't need the wfc files or .save dir for the qe calcs so we remove them at the end of each calc.
#		'''
		AFLOWpi.prep._addToAll(variedCalcs,'IMPORT','import shutil')
		AFLOWpi.prep._addToAll(variedCalcs,'IMPORT','import os')
		AFLOWpi.prep._addToAll(variedCalcs,'IMPORT','import glob')
		cleanupString="""number=1
while True:
    try:  
        wfcFile = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.wfc%d'%(ID,number))     
        os.remove(wfcFile)
        number+=1
    except Exception,e:
        break
try:
    shutil.rmtree(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.save'%ID))
except Exception,e:
    pass
"""
		AFLOWpi.prep._addToAll(variedCalcs,'CLEANUP',cleanupString)




############################################################################################################################
        



def _getCenter(origCalcs,fitVars=None,options=None,return_energy=False):
	minDict=AFLOWpi.pseudo._getMinimization(origCalcs,fitVars=fitVars,return_energy=return_energy) 
	returnDict=OrderedDict()
	for items in fitVars:
		returnDict[items]=minDict[0][items]
	if return_energy==True:
		try:
			minEnergy = minDict[0]['Energy']

			return returnDict,minEnergy
		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)
			return returnDict 
	else:
		return returnDict 

import __main__
# def _medE(calcs):
# 	calcs = __grabEnergyOut(calcs)
# 	energyList=[]
# 	for ID,oneCalc in calcs.iteritems():
# 		energyList.append(float(oneCalc['Energy']))
# 	return (sorted(energyList)[(len(energyList)/2)-1]+sorted(energyList)[(len(energyList)/2)])/2



# def _genetic(calcs):
# 	try:
# 		medE = __medE(calcs)

# 		survivors={}
# 		for ID,oneCalc in calcs.iteritems():
# 			if float(oneCalc['Energy'])>medE:
# 				survivors[ID]=oneCalc

		

# 		next_gen={}
# 		survivorsList = survivors.keys()
# 		while len(survivorsList)>0:
# 			choice1 = random.choice(survivorsList)
# 			survivorsList.remove(choice1)
# 			choice2 = random.choice(survivorsList)
# 			survivorsList.remove(choice2)
# 			pairs.append([choice1,choice2])

 #		for ID,oneCalc in survivors.iteritems()
				
#			if 
			
          # 	inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']
		# 	list1 = [{param:val} for param,val in inputDict.iteritems() if re.match(r'celldm',param)]
		# 	resultDict={}
		# 	for item in list1:
		# 		resultDict.update(item)


			
		# stepList=[steps for x in range(len(paramList))]
		# amountList=[amount for x in range(len(paramList))]
		
 		# variedCalcs=AFLOWpi.prep.varyCellParams(oneCalc,ID,param=paramList,amount=amountList,steps=stepList)
		# pass
#	except Exception,e:
#		AFLOWpi.run._fancy_error_log(e)


		




def _shiftGrid(calcs,outFile,fitVars=None,options=None,constraint=None,thresh=0.001,mult_jobs=True):
	'''sleep for anywhere from 0 to 10 seconds'''
	time.sleep (random.random()*10)
	'''if file semaphore exists, another process is working on submitting the grid so exit'''
	logging.debug('entering __shiftGrid')

	if os.path.exists(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE')):
		logging.debug('found GRID.SEMAPHORE. exiting __shiftGrid')
		sys.exit(0)
	
	else:
		'''if the file semaphore doesn't exist, place one and then, go ahead and shift grid'''
		with open(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE'),'w') as semaphoreFileObj:
			semaphoreFileObj.write('TRUE')


	try:

		'''
		in case we're 90% of the way through the walltime we don't want to risk the job possibly 
		getting killed half while it's submitting the jobs for the next iteration because it could 
		lead to double submission and a and a whole world of hurt
		'''



#		calcs[__main__.ID]['__walltime_dict__'][__main__.ID]['walltime']
		


		if AFLOWpi.prep._ConfigSectionMap('cluster','type').upper().strip() in ['PBS','UGE','SLURM']:
			'''
			in case we're 90% of the way through the walltime we don't want to risk the job possibly 
			getting killed half while it's submitting the jobs for the next iteration because it could 
			lead to double submission and a and a whole world of hurt
			'''



	#		calcs[__main__.ID]['__walltime_dict__'][__main__.ID]['walltime']
			
			if mult_jobs==True:
				walltime,startScript=AFLOWpi.run._grabWalltime(__main__.oneCalc,__main__.ID)
				if numpy.abs(time.time()-startScript)>walltime*0.90:
					'''in this case just resubmit this job'''
					AFLOWpi.run._submitJob(__main__.ID,__main__.oneCalc,__main__.__submitNodeName__,sajOverride=True)
					os.remove(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE'))		     
					sys.exit(0)			
                                        print 'quitting1'
					logging.debug('quitting1')
			else:

				'''in this case we need to find our way back to the calc that spawned the original 
				grid and restart that. it'll start back up and on it's first job it'll start up 
				where this script left off but will have the time to do it.'''

                                
				main_ID=os.path.basename(outFile).split('.')[0][1:]
				workdir=AFLOWpi.prep._ConfigSectionMap('prep','work_dir')
				mainOneCalc = AFLOWpi.prep._loadOneCalc(workdir,main_ID)
				startScript=mainOneCalc['__walltime_dict__']['start']
				walltime,startScript=AFLOWpi.run._grabWalltime(mainOneCalc,main_ID)

				if numpy.abs(time.time()-startScript)>walltime*0.90:
					logging.debug('quitting2')
					AFLOWpi.run._submitJob(main_ID,mainOneCalc,__main__.__submitNodeName__,sajOverride=True)

					os.remove(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE'))	      
					sys.exit(0)			


	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)


	try:

		try:
			if  type(list([4,5,6])) !=type(fitVars) or  type(tuple([4,5,6])) !=type(fitVars):
				fitVars=list(fitVars,)
		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)
		calcsEnergy = AFLOWpi.pseudo._grabEnergyOut(calcs)

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
					if constraints[1].upper() in ['A','B','C','ALPHA','BETA','GAMMA']:
						constraint_var_list.append('_AFLOWPI_PARAM_'+constraints[1].upper().strip()+'_')
					else:
						logging.warning('constraint: %s not valid. not including in minimization' % constraints)
			except Exception,e:
			    print e
			    constraint_type=None
			    constraint_var =None
		else:
			constraint_type=None
			constraint_var =None

		'''remove the constrained variables from the minimization'''
		for var in constraint_var_list:
			if var in fitVars:
				fitVars.remove(var)

		

		minEVal,minEnergy = AFLOWpi.pseudo._getCenter(calcs,fitVars=fitVars,options=options,return_energy=True)

		for ID,oneCalc in calcs.iteritems():
			if '__minEnergy__' not in oneCalc.keys():
				oneCalc['__minEnergy__']=[minEnergy]
			else:
				oneCalc['__minEnergy__'].append(minEnergy)

#		energy,varMatrix =  __getMatrices(calcs,fitVars=fitVars,options=options,Energy=False)		
		varMatrix =  AFLOWpi.pseudo._getMatrices(calcs,fitVars=fitVars,options=options,Energy=False)		

		gridDistanceDict=OrderedDict()

		for entry in range(len(fitVars)):
			gridDistanceDict[fitVars[entry]]=numpy.abs(varMatrix[0][entry]-varMatrix[1][entry])

		medianDict=OrderedDict()
		minDict=OrderedDict()
		maxDict=OrderedDict()

		for ID,oneCalc in calcs.iteritems():
			mainGridDir=os.path.dirname(os.path.dirname(oneCalc['_AFLOWPI_FOLDER_']))
			if '__gridMin_iteration__' not in oneCalc.keys():
				calcs[ID]['__gridMin_iteration__']=1
				iteration=1
			else:
				calcs[ID]['__gridMin_iteration__']=calcs[ID]['__gridMin_iteration__']+1
				iteration=calcs[ID]['__gridMin_iteration__']

			allVars = AFLOWpi.prep.extractvars(oneCalc['__refFile__'])
			refFile=oneCalc['__refFile__']
			break
		'''get min and max bounds for each variable in our combination space'''
		for item in range(len(fitVars)):
			'''for average the median is not in the center of the grid so we use mean'''
			medianDict[fitVars[item]]=numpy.median(numpy.asarray(varMatrix.T[item]))
			minDict[fitVars[item]]=numpy.min(varMatrix.T[item])
			maxDict[fitVars[item]]=numpy.max(varMatrix.T[item])
			
		'''check to see if value of minimum is within 1% of the bounds of each variable in our combination space'''
	
		close2MinList=[]
		close2MaxList=[]

		for fitVar,value in minEVal.iteritems(): 
			'''if we found the min to be within 2.5% an edge we need to reposition'''
			gFuncMin=numpy.abs(medianDict[fitVar]-minDict[fitVar])
			gFuncMax=numpy.abs(medianDict[fitVar]-maxDict[fitVar])
			close2MinList.append(numpy.abs(float(minEVal[fitVar])-float(minDict[fitVar]))<0.05*gFuncMin)
			close2MaxList.append(numpy.abs(float(minEVal[fitVar])-float(maxDict[fitVar]))<0.05*gFuncMax)


		'''
		now that we know if we're close edges we can shift the grid of test points centered around 1.95x
		the distance from the starting center point of the grid to where the minimization got stuck at
		'''
		shiftCenter=OrderedDict()
		logging.debug('close2MinList %s ' % close2MinList)
		logging.debug('close2MaxList %s ' % close2MaxList)


		for item in range(len(fitVars)):
			if close2MinList[item]:
				'''if it's at the min bound shift should make the center go down'''
				shiftCenter[fitVars[item]] = -numpy.abs(medianDict[fitVars[item]]-minDict[fitVars[item]])
			elif close2MaxList[item]:
				'''if it's at the max bound shift should make the center go up'''
				shiftCenter[fitVars[item]] =  numpy.abs(medianDict[fitVars[item]]-maxDict[fitVars[item]])
			else:
				shiftCenter[fitVars[item]]=0


		'''
		identify the variables assiciated with the fitVars by looking at the reference file
		(this isn't the best way to do this but works for the case of celldm)
		if needed it'll be made more robust in the future
		'''
		assocInputVars=OrderedDict()
		try:
			refFileDict = AFLOWpi.retr._splitInput(refFile)
		except:
			pass
		
		for k,v in refFileDict.iteritems():
			for l,m in v.iteritems():
				if m in fitVars:
					assocInputVars[l]=m

		'''now we can generate the shifted set from the old set'''
		shiftedCalcs=copy.deepcopy(calcs)

		'''check to see if our old min is close to the new min'''
		quit_bool=True
		centerShiftDict=OrderedDict()
		for item in fitVars:
			centerShift=medianDict[item]-minEVal[item]
			centerShiftDict[item]=centerShift
			if abs(centerShift)>thresh:
				quit_bool=False
		'''write out the hessian for bulk modulus if we're done minimizing'''
		if quit_bool==True:
			try:
				AFLOWpi.pseudo._getMin(calcs,fitVars=fitVars,bulk_modulus=True)
			except Exception,e:
				AFLOWpi.run._fancy_error_log(e)
				

		shiftType='shift'
		if True not in close2MaxList and True not in close2MinList:
			if quit_bool==False:
				shiftType='shrink'
			else:
				shiftType='finish'

		with open(os.path.join(os.path.dirname(outFile),'GridEnergyLog.log'),'a+') as myfile:
			myfile.write('\n{:=^31}\n'.format(''))
			myfile.write('{:=^31}\n'.format('iteration %02d'%iteration))
			myfile.write('{:=^31}\n'.format(''))
#			myfile.write('----------------------------------\n')
			myfile.write('{:<30}'.format('Status |   %s'%(shiftType.upper())))
			myfile.write('|\n')

			myfile.write('Energy |')
                        try:
                            myfile.write('{: 19.12f}'.format(float(minEnergy)))
                        except:
                            print minEnergy
			myfile.write(' Ry|\n')


#			myfile.write('\n')


			paramsPrint=[x.split('_')[-2].upper() for x in minEVal.keys()]
			myfile.write('--------')
			myfile.write('-----------------------\n'*len(paramsPrint))
			myfile.write('Param  | ')

			for item in range(len(paramsPrint)):
				if paramsPrint[item] in ['A','B','C']:
					paramsPrint[item]+=' (Bohr)'

				if paramsPrint[item] in ['ALPHA','BETA','GAMMA']:
					paramsPrint[item]+=' (Deg.)'
				myfile.write('{0:*^18}   |'.format(paramsPrint[item]))
			myfile.write('\n')
			myfile.write('--------')
			myfile.write('-----------------------'*len(paramsPrint))
			myfile.write('\n')
			myfile.write('Median |  ')
			for k,v in minEVal.iteritems():
				myfile.write(' {: 16.11f}   |'.format(numpy.around(medianDict[k],decimals=10)))
			myfile.write('\n')
			myfile.write('Min Val|  ')
			for k,v in minEVal.iteritems():
				myfile.write(' {: 16.11f}   |'.format(numpy.around(minEVal[k],decimals=10)))
			myfile.write('\n')
			myfile.write('Delta  |  ')
			for k,v in minEVal.iteritems():
				myfile.write(' {: 16.11f}   |'.format(numpy.around(-1.0*centerShiftDict[k],decimals=10)))
			myfile.write('\n')



#		calcs=AFLOWpi.retr.grabEnergyOut(calcs)

		fileNameStr=os.path.join(mainGridDir,'energyLandscape_iteration_%02d.pdf' % iteration)
		if not os.path.exists(fileNameStr):
			filenameStr='./energyLandscape_iteration_%02d.pdf' % iteration
		try:
			if len(fitVars)==2:
				xTitle=fitVars[0].split('_')[-2]
				xTitle+=' (Bohr)'
				yTitle=fitVars[1].split('_')[-2]
				yTitle+=' (Bohr)'
				AFLOWpi.plot.interpolatePlot(calcs,fitVars[0],fitVars[1],xaxisTitle=xTitle,yaxisTitle=yTitle,fileName=fileNameStr,circle_min=True,text_min=True,vhline_min=True,delta_min=False)		

				AFLOWpi.plot.interpolatePlot1D(calcs,fitVars[0],xaxisTitle=xTitle,fileName='X_'+fileNameStr,circle_min=True) 
				AFLOWpi.plot.interpolatePlot1D(calcs,fitVars[1],xaxisTitle=xTitle,fileName='Y_'+fileNameStr,circle_min=True) 
                                raise SystemExit
			if len(fitVars)==1:
				xTitle=fitVars[0].split('_')[-2]
				AFLOWpi.plot.interpolatePlot1D(calcs,fitVars[0],xaxisTitle=xTitle,fileName=fileNameStr,circle_min=True) 

		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)
                
		'''if we find the min in the bounds this is how much we'll shrink the grid by'''
		shrinkVarDict={}
		for item in fitVars:
			try:
				shrinkVarDict[item]=1.0-numpy.abs(centerShiftDict[item]/(medianDict[item]-maxDict[item]))
				if shrinkVarDict[item]>0.75:
					'''put a limit so if shift is very small
					we don't scale down 1000x or something'''
					shrinkVarDict[item]=0.75

			except Exception,e:
				AFLOWpi.run._fancy_error_log(e)
				shrinkVarDict[item]=0.49


                ibrav=0
                try:
                    for ID,oneCalc in shiftedCalcs.iteritems():
                        inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
                        ibrav=int(inputDict['&system']['ibrav'])
                        break
                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)

		'''if the old min is close to the new min we won't resubmit and instead   '''
		'''and instead we will look copy the results back to the parent directory '''


		for ID,oneCalc in shiftedCalcs.iteritems():
			for item in fitVars:

				'''if it's not close to any bounds we want to generate 
				an input file with the optimal coordinates'''



				if True not in close2MaxList and True not in close2MinList:
					logging.info('minimum found within grid bounds. shifting %s by %s and shrinking grid'%(item,centerShiftDict[item]))
				else:
					'''if we don't find it in the bounds we've already shifted so stay with the
					   same size grid just shifted'''
					logging.info('minimum not found within grid bounds. shifting %s by %s'%(item,centerShiftDict[item]))
					
				'''just set all of the calcs to the min (we won't end up saving
				   these inputs. because of quit_bool the calcs will make an output
				   file in the parent directory with the params of the first calc in the
				   loop (which will then have the min params in it)'''
				if not quit_bool:
					'''
					shrink grid by shifting all points closer to new center scale the
					shrinking by how much we've shifted compared to grid width
					'''
					shiftedCalcs[ID][item]=shiftedCalcs[ID][item]-centerShiftDict[item]
					shrinkPrint=1.0-shrinkVarDict[item]
					logging.debug('shrink factor for %s=%s'%(item,shrinkPrint))

					centDist=(shiftedCalcs[ID][item]-minEVal[item])
					scaledCentDist=centDist*shrinkVarDict[item]
					logging.debug('before: %s' %shiftedCalcs[ID][item])
					logging.debug('center: %s' %minEVal[item])
					shiftedCalcs[ID][item]-=scaledCentDist
					'''round the numbers to 9 decimal places so we can reform for plotting'''
					logging.debug('after: %s' %shiftedCalcs[ID][item])


				else:
					shiftedCalcs[ID][item]=minEVal[item]				       

						

			if constraint!=None:
				for constr in range(len(constraint_type_list)):
					if constraint_type_list[constr]=='volume':
						"""NEED TO BE GENERALIZED NEED TO BE GENERALIZED NEED TO BE GENERALIZED"""
						vol = float(AFLOWpi.retr.getCellVolume(oneCalc,ID))

						"""NEED TO BE GENERALIZED NEED TO BE GENERALIZED NEED TO BE GENERALIZED"""
						constr_var = constraint_var_list[constr]

						volDiv=1
						for param in ['A','B','C']:
							try:
								param_name = '_AFLOWPI_PARAM_%s_'%param
								if param_name!=constr_var:
									volDiv*=oneCalc[param_name] 
							except Exception,e:
								AFLOWpi.run._fancy_error_log(e)

						for param in ['ALPHA','BETA','GAMMA',]:
							if param in oneCalc.keys():
								param_name = '_AFLOWPI_PARAM_%s_'%param
								if param_name!=constr_var:
									temp_angle_param=[param_name] 
									volDiv*=numpy.sin(temp_angle_param*(numpy.pi/180.0))
						remainder=vol/volDiv

						if constr_var.upper() in ['_AFLOWPI_PARAM_ALPHA_','_AFLOWPI_PARAM_BETA_','_AFLOWPI_PARAM_GAMMA_',]:
							oneCalc[constr_var.upper().strip()]=numpy.arcsin(remainder)*(180.0/numpy.pi)
						if constr_var.upper() in ['_AFLOWPI_PARAM_A_','_AFLOWPI_PARAM_B','_AFLOWPI_PARAM_C_',]:
							oneCalc[constr_var.upper().strip()]=remainder

#                                        if constraint_type_list[constr]=='fixed':
#                                            oneCalc[constr_var.upper().strip()]=oneCalc[constr_var.upper().strip()]



		for ID,oneCalc in shiftedCalcs.iteritems():

			
			inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
			inputDict['&control']['restart_mode']="'from_scratch'"
						

			'''convert A,B,C,alpha,beta,gamma to celldm'''
			for l,m in inputDict['&system'].iteritems():

					if l=='celldm(1)':
						inputDict['&system']['celldm(1)']=oneCalc['_AFLOWPI_PARAM_A_']
                                        try:
                                            if l=='celldm(2)':
						inputDict['&system']['celldm(2)']=oneCalc['_AFLOWPI_PARAM_B_']/oneCalc['_AFLOWPI_PARAM_A_']
                                        except:
                                            pass
                                        try:
                                            if l=='celldm(3)':
						inputDict['&system']['celldm(3)']=oneCalc['_AFLOWPI_PARAM_C_']/oneCalc['_AFLOWPI_PARAM_A_']
                                        except:
                                            pass
                                        try:
                                            if l=='celldm(4)':
						newcelldm=numpy.arccos(oneCalc['_AFLOWPI_PARAM_ALPHA_']*180.0/numpy.pi)
						inputDict['&system']['celldm(4)'] = newcelldm
                                        except:
                                            pass
                                        try:
                                            if l=='celldm(5)':
						newcelldm=numpy.arccos(oneCalc['_AFLOWPI_PARAM_BETA_']*180.0/numpy.pi)
						inputDict['&system']['celldm(5)'] = newcelldm
                                        except:
                                            pass
                                        try:
                                            if l=='celldm(6)':
						newcelldm=numpy.arccos(oneCalc['_AFLOWPI_PARAM_GAMMA_']*180.0/numpy.pi)
						inputDict['&system']['celldm(6)'] = newcelldm
                                        except:
                                            pass

			if quit_bool==False:# in close2MaxList or True in close2MinList:
				
				newInputStr = AFLOWpi.retr._joinInput(inputDict)
				shiftedCalcs[ID]['_AFLOWPI_INPUT_']=newInputStr


				'''set the status of this calc set to not started or completed'''
				shiftedCalcs[ID]['__status__']=OrderedDict()
				shiftedCalcs[ID]['__status__']["Start"]=False
				shiftedCalcs[ID]['__status__']["Complete"]=False
				shiftedCalcs[ID]['__status__']["Restart"]=0
				shiftedCalcs[ID]['__status__']["Error"]="None"
				'''set the execCounter to 0 so when these calcs get resubmitted they'll run fresh'''
				shiftedCalcs[ID]['__execCounter__']=0
#				for var in range(len(fitVars)):
#					shiftedCalcs[ID][fitVars[var]]=float(shiftedCalcs[ID][fitVars[var]])+float(shiftCenter[fitVars[var]])
				AFLOWpi.prep._saveOneCalc(shiftedCalcs[ID],ID)



				with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % ID),'w') as outInfileObj:
					outInfileObj.write(newInputStr)
				'''since we shifted the vals in the input file we have to make sure we shift the vals in oneCalc'''
				

				
			else:
				'''
				if we found the minimum was not on the border of the grid we save
				a copy of the input with the minimized parameters in the AFLOWpi folder
				and then quit because we don't have to submit anymore jobs.
				'''

#				cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(newInputStr)
				'''we need alat to scale matrix down like it is in espresso output'''
#				alat = float(inputDict['&system']['celldm(1)'])
#				cellParamMatrix/=alat
#				matrixString = AFLOWpi.retr._cellMatrixToString(cellParamMatrix)
#				inputDict['CELL_PARAMETERS']=OrderedDict()
#				inputDict['CELL_PARAMETERS']['__content__']=matrixString
				'''format like espresso output'''
#				inputDict['CELL_PARAMETERS']['__modifier__']='(alat = %f)' % alat
#				atomPOS = inputDict['ATOMIC_POSITIONS']['__content__']
#				detachedPos,flags=AFLOWpi.retr.detachPosFlags(atomPOS)
#				inputDict['ATOMIC_POSITIONS']['__content__']=detachedPos
				newInputStr = AFLOWpi.retr._joinInput(inputDict)
                                inputDict = AFLOWpi.retr._splitInput(newInputStr)
                                
                                outFile_split=list(os.path.split(outFile))
                                dest_ID=outFile_split[-1].split('.')[0][1:]

                                dest_calc = AFLOWpi.prep._loadOneCalc(outFile_split[0],dest_ID)

                                dest_calc_dict=AFLOWpi.retr._splitInput(dest_calc['_AFLOWPI_INPUT_'])
                                for k,v in inputDict['&system'].iteritems():
                                    dest_calc_dict['&system'][k]=v

				newInputStr = AFLOWpi.retr._joinInput(dest_calc_dict)
                                new_in=os.path.join(os.path.dirname(outFile),"%s.in"%dest_ID)
				with open(new_in,'w') as outInfileObj:
					outInfileObj.write(newInputStr)
                                

                                dest_calc['_AFLOWPI_INPUT_']=newInputStr
                                AFLOWpi.prep._saveOneCalc(dest_calc,dest_ID)

				'''remove file semaphore'''
				try:
					os.remove(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE'))
				except Exception,e:
					AFLOWpi.run._fancy_error_log(e)

				
				if mult_jobs==True:
					return True
				else:
					oldConfig = os.path.join(os.path.dirname(os.path.dirname(outFile)),'AFLOWpi','CONFIG.config')
					AFLOWpi.prep._forceGlobalConfigFile(oldConfig)
					return True

				logging.debug('exiting __shiftGrid')
				


		'''if we didn't find the minimum we submit the shifted grid to run'''
                try:
                    sub_node_name = __main__.__submitNodeName__
                except:
                    sub_node_name = socket.gethostname()
		globals()['__submitNodeName__'] = sub_node_name

		'''remove file semaphore'''
		try:
			os.remove(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE'))
		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)



		
		
		'''clear out all wfc and .save scratch fromt he dirs'''
		for ID,oneCalc in reversed(shiftedCalcs.items()):
			number=1
			while True:
				try:  
					wfcFile = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.wfc%d'%(ID,number))     
					os.remove(wfcFile)
					number+=1
				except Exception,e:
					break
			try:
				shutil.rmtree(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.save'%ID))
			except Exception,e:
				pass



		'''if we are submitting the grid calc jobs separately or one big job'''
		if mult_jobs==True:
                    oneJobBool=False
                    sajOver=True
		else:
                    sajOver=False
                    oneJobBool=True
		'''submit in reverse order because calcs later in the orderedDict are more likely'''
		'''to be larger cells than those at the beginning'''

		for ID_new,oneCalc_new in reversed(shiftedCalcs.items()):
			AFLOWpi.run._submitJob(ID_new,oneCalc_new,sub_node_name,forceOneJob=oneJobBool,sajOverride=sajOver)
			# try:
			# 	fileList = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'./_*wfc*')

			# 	logging.debug('FILELIST %s' %fileList)
			# 	for fileName in fileList:
			# 		try:
			# 			os.system('rm  %s' %fileName )


			# 		except Exception,e:
			# 			AFLOWpi.run._fancy_error_log(e)

			# 	shutil.rmtree(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.save'%ID))

			# except Exception,e:
			# 	AFLOWpi.run._fancy_error_log(e)

		logging.debug('exiting __shiftGrid')
		return False




	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
		try:
			os.remove(os.path.join(os.path.dirname(outFile),'GRID.SEMAPHORE'))
			raise SystemExit
		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)
			raise SystemExit

		


def _getMinimization(origCalcs,fitVars=None,options=None,return_energy=False,minimize_var="Energy"):
    '''
    Looks through the dictionary of dictionaries that you supply it and finds the minimum energy
    for the different cutoff/kpoint choices. outputs a list of dictionaries with the information
    about the cutoffs and the value of the variables that gives you a minimum energy.

    Arguments:
          origCalcs (list): a list of Dictionary of Dictionaries of the calculations.

    Keyword Arguments:
          fitVars (list): tuple of variables that you want to minimize energy with respect to.
          options (dict): additional options passed to L-BFGS-B minimization 
                          (see scipy documentation for scipy.optimize.minimize for more details)
          return_energy (bool): If true add the min energy to the resultList list for each
          minimize_var (str): which key to minimize on in oneCalc

    Returns:
          resultList (list): list if minumum energies and the parameters of the set that they 
                             correspond to.
                         
    '''
    try:
	    key = AFLOWpi.pseudo._getCutOffs(origCalcs)
    except Exception,e:
            AFLOWpi.run._fancy_error_log(e)
	    key=['']
    manyCalcs=AFLOWpi.pseudo._splitCalcs(origCalcs,key)
    
#    print "Entering getMinimization"
    if fitVars==None:
        print 'ERROR: You need to specify variables for the fit and minimization'
	logging.warning('ERROR: You need to specify variables for the fit and minimization')
        return
    
    resultList = []
    for calcs in manyCalcs:

        '''
        actually grabs your energy from each of the calculations
        and puts it into a dictionary of dictionaries containing 
        all the parameters we will need to calculate the minimum
        '''
        energyDict = AFLOWpi.pseudo._grabEnergyOut(calcs)


	for ID,oneCalc in energyDict.iteritems():
		sampleDict=copy.deepcopy(oneCalc)
		break
        resultDict= OrderedDict()       

	'''
	checks to see if you are scaling ecutrho to ecutwfc via the '_AFLOWPI_DUAL_' variable
	'''

	for cutoff in key:
		if cutoff=='_AFLOWPI_ECUTR_':
			if sampleDict.has_key('_AFLOWPI_DUAL_'):
				resultDict['_AFLOWPI_DUAL_']=sampleDict['_AFLOWPI_DUAL_']
		else:
			resultDict[cutoff]=sampleDict[cutoff]

	energyMatrix,varMatrices = AFLOWpi.pseudo._getMatrices(energyDict,fitVars=fitVars,options=options)
	try:
		minEnergy,minValue = AFLOWpi.pseudo._getMin(energyDict,fitVars=fitVars,options=options,minimize_var=minimize_var)
		minValue = OrderedDict(resultDict.items()+ minValue.items())
		
		if return_energy==True:
			minValue['Energy']=minEnergy
		resultList.append(minValue) 
	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
		pass

    
#    print "Exiting getMinimization"
    return resultList
            

from operator import itemgetter

def _plotOne(plots,labs,fig,entry,key,xaxis,pltTitle=None,rename=None,entryNum=0,maxs=None,mins=None):
	'''
	takes in a list of dictionaries of dictionaries of calculations and generates
	a plot with the x axis being some value in the list 'key' and splits the calculations
	and plots them with each plot being a unique combination of the items in key that are not	
'xaxis'
	Arguments:
	      entry (list): list of dictionaries of dictionaries of calculations
	      key (list): a list of cutoff variables used in the calculations
   	      xaxis (str): the cutoff that you choose to be the x axis in your plots

        Keyword Arguments:
	      plotTitle -- title of the plots (default: None)
	'''

	"""
	takes a dictionary which is called renamed whose keys are the current variable names 
	and the values are what you want to replace them with. this allows to customize axis 
	labels on the plot
	"""
	
	filename='PSEUDOTEST'+str(entryNum)
#	for cutoff in range(len(key)):
#		if key[cutoff]=='_AFLOWPI_KPOINTS_':
#			filename+='_%s' % ('-'.join(entry[0][key[cutoff]].split()[:3]))
#		else:
#			filename+='_%s' % str(entry[0][key[cutoff]])
	lineName = ''
	for cutoff in key:
		if cutoff in entry[0].keys():
			lineName+=' '+cutoff+' '+str(entry[0][cutoff])
#	print lineName

	if pltTitle==None:
		pltTitle=filename
	
	if rename!=None:
		for oneCalc in entry:
			entryCopy=copy.deepcopy(oneCalc)
			for ID,value in entryCopy.iteritems():
				if ID in rename.keys():
					oneCalc[rename[ID]]=value
					del oneCalc[ID]
		for ID in range(len(key)):
			if key[ID] in rename.keys():
				
				key[ID]=rename[key[ID]]
				
		if xaxis in rename.keys():
			xaxis=rename[xaxis]
				
		
        for oneSet in entry:
		fitVarList = []
		for k,v in oneSet.iteritems():
			if k not in key:
				if k != xaxis:
					fitVarList.append(k)
        
	
                
        color_cycle=[ 'm','c','r', 'g', 'b', 'y',]
	minMaxList=[]
        for var in range(len(fitVarList)):
            for value in entry:  
		  minMaxList.append(value[fitVarList[var]])

	maxY = max(minMaxList)
	minY = min(minMaxList)
	rangeY=maxY-minY

	maxY+=1.1*rangeY
	minY-=1.1*rangeY


        for var in range(len(fitVarList)):

            varYaxis=[]
            varXaxis=[]
            for value in entry:  
                varXaxis.append(float(value[xaxis]))
                varYaxis.append(value[fitVarList[var]])

	    ax = plt.subplot(len(fitVarList),1,var+1)
 
	    lns = plt.plot(varXaxis,varYaxis,marker='.',linestyle='-',color=color_cycle[entryNum%len(color_cycle)],label=lineName)
	    if var==0:
		    plt.title(r'%s' % pltTitle)

            rangeY=numpy.abs(maxs[var]-mins[var])
            if rangeY < 0.50:
                rangeY=0.1
            if maxs[var]<0:
                maxs[var]-=rangeY*-0.49
            else:
                maxs[var]-=rangeY*-0.49
            
            if mins[var]<0:
                mins[var]-=rangeY*-0.49
            else:
                mins[var]-=rangeY*0.49
            

	    plt.ylim(mins[var],maxs[var])
	    plots.append(lns)
	    ax.set_ylabel(fitVarList[var])                
	    if var==len(fitVarList)-1:
		    ax.set_xlabel(xaxis)		

#	    labs = [l.get_label() for l in lns]
	    labs.append(lineName)
#	    print plots
	    if var==0:
		    plt.legend()
	
	boxColors = ['darkkhaki','royalblue']

	try:
		key.remove(xaxis)
	except ValueError,e:
		pass

       
	resDir = './'
#	if not os.path.exists(resDir):
#		os.mkdir(resDir)
	
#	plt.savefig(os.path.join(resDir,'%s.pdf' % filename),bbox_inches='tight')
	return labs,plots
	    
def plot(resultList,xaxis='',xtitle=None,ytitle=None,title=None,rename=None,file_name=''):
    '''
    takes in a list of dictionaries of dictionaries of calculations and generates
    a plot with the x axis being some value in the list 'key' and splits the calculations
    and plots them with each plot being a unique combination of the items in key that are not
    'xaxis'

    Arguments:
          resultList (list): list of dictionaries of dictionaries of calculations
          xaxis (str): the keyword in oneCalc that you choose to be the x axis in your plots

    Keyword Arguments:
          xtitle (str): title of the x axis of the plot (default: None)
          ytitle (str): title of the y axis of the plot (default: None)
          plotTitle (str): title of the plots (default: None)
          rename (dict): a mapping of the names of the keywords of whose
                         values used to generate the plot to some other name.
                         ex. {' _AFLOWPI_ECUTW_':'wavefunction cutoff'}
          file_name (str): use this instead of "PT_RESULTS.pdf" as filename of plot

    Returns:
          None
    '''


    print "Entering generatePlot"
    resultListDict=OrderedDict()
    [resultListDict.update({str(i):resultList[i]}) for i in range(len(resultList))]

    key = AFLOWpi.pseudo._getCutOffs(resultListDict)
    
    if resultList[0].has_key('_AFLOWPI_DUAL_'):
	    key.append('_AFLOWPI_DUAL_')
    
    value=[]
    splittingForPlot = list([entry for entry in key if entry != xaxis])    
    
    splitResults =  AFLOWpi.pseudo._splitCalcs(resultListDict,splitVars=splittingForPlot)
#    print splittingForPlot
#    print splittingForPlot
#    print splittingForPlot
    calcsFix = []
    '''
    converts the dictionary of dictionaries back into a list of dictionaries
    '''
    for entry in splitResults:

        for i in entry.keys():
            entry[i][xaxis]=int(entry[i][xaxis])

        calcsFix.append(entry.values())
    
    '''
    plots the list of dictionaries that have been split on the cutoffs
    '''
    entryNum=0
    width = 20
    height = 14
    fig = pylab.figure(figsize=(width, height))
    labs=[]
    plots=[]
    color_cycle=['r', 'g', 'b', 'y','c', 'm']

    lim_dict=OrderedDict()
    for entry in calcsFix:
        for k in entry:
            for j in k.keys():
                if j!=xaxis and j not in splittingForPlot:
                    try:
                        lim_dict[j].append(float(k[j]))
                    except:
                        lim_dict[j]=[]
    max_list=[]
    min_list=[]
    for k,v in lim_dict.iteritems():
        
        max_list.append(max(v))
        min_list.append(min(v))

    for entry in calcsFix:

        entrySorted = sorted(entry, key=lambda x: [x[cutoff] for cutoff in key])


#        if title==None:

        labs,plots = AFLOWpi.pseudo._plotOne(plots,labs,fig,entrySorted,splittingForPlot,xaxis,pltTitle=title,rename=rename,entryNum=entryNum,maxs=max_list,mins=min_list)

        
#	print labs
#	plt.figtext(0.10, 0.98-0.035*entryNum, labs[entryNum], 
#		    backgroundcolor=color_cycle[((entryNum+1)%len(color_cycle))-1], color='black', weight='roman',
#		    size='large')
        entryNum+=1
    if file_name!='':
        filename=file_name
    else:
        filename = 'PT_Results'
    
    resDir = './'
    plt.savefig(os.path.join(resDir,'%s.pdf' % filename),bbox_inches='tight')
    
    print "Exiting generatePlot"
