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

import logging
import AFLOWpi
import os
import __main__
import numpy
import time
import configparser
import glob 
import sys 
import re


def _add_subset_to_daemon_log_list(addition_list,log_name):
    if  AFLOWpi.prep._ConfigSectionMap('cluster','daemon').lower()=='true':
        AFLOWpi.run._submit_log_append(addition_list,log_name)


def _one_test_build(oneCalc,ID,build_command,subset_name='SUBSET',merge_oneCalc=True,keep_name=False,config=None,clean_input=True):
    if config==None:
        config=oneCalc['_AFLOWPI_CONFIG_']

    intoInit={'PROJECT':subset_name,'SET':'','workdir':oneCalc['_AFLOWPI_FOLDER_'],'config':config}
    fake_session_keys = AFLOWpi.prep.init(**intoInit)

    input_strings=None
    _locals=locals()
    exec('input_strings=%s'%build_command,globals(),_locals)
    input_strings = _locals["input_strings"]

    if merge_oneCalc==True:
        varied_calcs = AFLOWpi.prep.calcFromFile(fake_session_keys,input_strings,reffile=oneCalc['_AFLOWPI_INPUT_'],workdir=oneCalc['_AFLOWPI_FOLDER_'],keep_name=keep_name,clean_input=clean_input)
    else:
        varied_calcs = AFLOWpi.prep.calcFromFile(fake_session_keys,input_strings,workdir=oneCalc['_AFLOWPI_FOLDER_'],keep_name=keep_name,clean_input=clean_input)
    

    return varied_calcs
###############################################################################################################

###############################################################################################################

def prep_split_step(calcs,subset_creator,subset_tasks=[],mult_jobs=False,substep_name='SUBSET',keep_file_names=False,clean_input=True,check_function=None,fault_tolerant=False,pass_through=False):


        
#####################################################################
        AFLOWpi.run._skeletonRun(calcs) 
        if check_function!=None:
            check_function=repr(check_function)

        for ID,oneCalc in list(calcs.items()):

                oneCalc['__splitCounter__']=0

                execString='''if oneCalc['__execCounter__']<=%s:
    ''' % oneCalc['__execCounterBkgrd__']
                execString+='''
     oneCalc,ID = AFLOWpi.prep.construct_and_run(__submitNodeName__,oneCalc,ID,build_command="""%s""",subset_tasks=%s,mult_jobs=%s,subset_name='%s_%s',keep_file_names=%s,clean_input=%s,check_function=%s,fault_tolerant=%s,pass_through=%s)

''' % (subset_creator,repr(subset_tasks),mult_jobs,ID,substep_name,keep_file_names,clean_input,check_function,fault_tolerant,pass_through)  
#''' % (subset_creator,exit_command,repr(subset_tasks),mult_jobs)  
                oneCalc['__execCounterBkgrd__']+=1
                AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN', execString)

                execString='''if oneCalc['__execCounter__']<=%s:
     ''' % oneCalc['__execCounterBkgrd__']
                execString+='''del oneCalc["__CRAWL_CHECK__"]
     oneCalc["__execCounter__"]+=1
     AFLOWpi.prep._saveOneCalc(oneCalc,ID)
'''
                oneCalc['__execCounterBkgrd__']+=1
                AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN', execString)
#####################################################################                

        
        return calcs

#####################################################################################################################




########################################################################################################################################################################################################################################

def construct_and_run(__submitNodeName__,oneCalc,ID,build_command='',subset_tasks=[],fault_tolerant=False,mult_jobs=True,subset_name='SUBSET',keep_file_names=False,clean_input=True,check_function=None,pass_through=False):

#        subset_name = ID + "_" + subset_name

        sub_path=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name)
        if not os.path.exists(sub_path):
            os.mkdir(sub_path)
                              
        '''this is a check to see if we're restarting when mult_jobs==True'''
        checkBool=False
        if '__CRAWL_CHECK__' in list(oneCalc.keys()):
            if oneCalc['__CRAWL_CHECK__']==ID:
                checkBool=True


        chain_index=1
        try:

            chain_index=oneCalc['__chain_index__']
            #           AFLOWpi.prep._passGlobalVar('__TEMP__INDEX__COUNTER__',oneCalc['__TEMP__INDEX__COUNTER__'])
            AFLOWpi.prep._passGlobalVar('__TEMP__INDEX__COUNTER__',chain_index)
        except Exception as e:
                AFLOWpi.run._fancy_error_log(e)
        
        chain_logname='step_%02d'%1
        logging.debug(chain_logname)
        logging.debug('CHECKBOOL:%s'%checkBool)
        if checkBool==False:

            AFLOWpi.prep._from_local_scratch(oneCalc,ID)

            if check_function==None:
                complete_function='True'
            else:
                complete_function=check_function
            #block this prep from being run again.
            oneCalc['__CRAWL_CHECK__']=ID
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)

            outFile=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in')
            command = '''
         completeBool=%s
         if completeBool:
            workdir = '../../'
            mainOneCalc = AFLOWpi.prep._loadOneCalc(workdir,'%s')
            AFLOWpi.prep._swap_walltime_logs('%s',mainOneCalc,oneCalc,ID)
            AFLOWpi.run._submitJob('%s',mainOneCalc,__submitNodeName__,forceOneJob=True)

''' % (complete_function,ID,ID,ID)
########y########################################################################################################
            try:
                os.mkdir(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi'))
            except:
                pass

            newConfigPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi','CONFIG.config')
            config = configparser.RawConfigParser()
            config.read(oneCalc['_AFLOWPI_CONFIG_'])
            config.set('prep', 'work_dir', oneCalc['_AFLOWPI_FOLDER_']) 

            if config.has_section('cluster'):
                if config.has_option('cluster','job_template'):
                    try:

                        qsub_temp_ref = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi','CLUSTER.ref')
                        qsubSub='''cd .*%s\npython .*%s''' % (os.path.basename(oneCalc['_AFLOWPI_FOLDER_']),os.path.join(os.path.basename(oneCalc['_AFLOWPI_FOLDER_']),'_'+ID+'.py'))

                        qsubSub_reg = re.compile(qsubSub)

                        with open(oneCalc['__qsubFileName__'],'r') as qsub_pre_trans:
                            qsub_string = qsub_pre_trans.read()

                        qsub_string = qsubSub_reg.sub('',qsub_string)

                        with open(qsub_temp_ref,'w') as qsub_post_trans:
                            qsub_post_trans.write(qsub_string)

                        config.set('cluster', 'job_template',qsub_temp_ref) 
                    except Exception as e:
                        AFLOWpi.run._fancy_error_log(e)

            with open(newConfigPath,'w') as fileWrite:    
                config.write(fileWrite)

################################################################################################################
            calc_subset = AFLOWpi.prep._one_test_build(oneCalc,ID,build_command,subset_name=subset_name,keep_name=keep_file_names,config=newConfigPath,clean_input=clean_input)

            AFLOWpi.prep.runAfterAllDone(calc_subset,command,faultTolerant=fault_tolerant)


            '''if we are submitting the grid calc jobs separately or one big job'''



            for task in subset_tasks:                
                _locals=locals()
                exec(task,globals(),_locals)


            for ID_sub,oneCalc_sub in list(calc_subset.items()):
                set_complete_string='''oneCalc['__status__']['Complete']=False
AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''
                AFLOWpi.prep._addToBlock(oneCalc_sub,ID_sub,'LOCK',set_complete_string) 

                set_complete_string='''oneCalc['__status__']['Complete']=True
AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''
                AFLOWpi.prep._addToBlock(oneCalc_sub,ID_sub,'SUBMITNEXT',set_complete_string) 

            '''submit in reverse order because calcs later in the orderedDict are more likely'''
            '''to be larger cells than those at the beginning'''                
 #           invert_bool=True
            '''if we're almost at the end of the walltime don't try to submit'''
            walltime,startScript=AFLOWpi.run._grabWalltime(oneCalc,ID)
            try:
                bn = os.path.basename(oneCalc['_AFLOWPI_FOLDER_'])
                subset_logs='../../%s/%s/AFLOWpi/calclogs/%s.log'%(bn,subset_name,chain_logname)
                AFLOWpi.prep._add_subset_to_daemon_log_list([subset_logs],'../AFLOWpi/submission_daemon/log_list.log')
            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)
                

            #keep track of time in main script as it loops in case
            #we are running serial jobs and the walltime runs out
        else:
            subset_config = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],subset_name,'AFLOWpi','CONFIG.config')
            calc_subset=AFLOWpi.prep.loadlogs(subset_name,'',chain_logname,config=subset_config)
            logging.debug(list(calc_subset.keys()))
            try:
                walltime,startScript=AFLOWpi.run._grabWalltime(oneCalc,ID)
            except:
                pass
            #exit if all calcs are done (needed for local mode)
            try:
                for k,v in list(calc_subset.items()):
                    oneCalc_sub=v
                    ID_sub=k
                    break
                if AFLOWpi.prep._checkSuccessCompletion(oneCalc_sub,ID_sub,faultTolerant=fault_tolerant):
                    return oneCalc,ID
            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

#################################################################################################
        for ID_new,oneCalc_new in list(calc_subset.items()):
            try:
                calc_subset[ID_new]['__walltime_dict__']=oneCalc['__walltime_dict__']
                AFLOWpi.prep._saveOneCalc(oneCalc_new,ID_new)
            except Exception as e:
                print(e)
                pass
#################################################################################################            
        if mult_jobs==True:
            oneJobBool=False
            sajO=True
            
        else:
            oneJobBool=True
            sajO=False
            AFLOWpi.prep._return_to_main_pipeline(calc_subset,oneCalc,ID)

        try:
                last=len(calc_subset)
                for ID_new,oneCalc_new in list(calc_subset.items()):

                    last-=1 
                    if last==0:
                        #to make sure this doesn't try to run again
                        oneCalc['__execCounter__']+=1
                        oneCalc['prev'].append(ID)
                        AFLOWpi.prep._saveOneCalc(oneCalc,ID)

                        oneJobBool=True
                        sajO=False

                    AFLOWpi.run._submitJob(ID_new,oneCalc_new,__submitNodeName__,forceOneJob=oneJobBool,sajOverride=sajO)
                if not pass_through:
                    sys.exit(0)                    

        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)
            sys.exit(0)
            
        
            
        return oneCalc,ID


def _return_to_main_pipeline(calc_subset,oneCalc,ID):
    '''
    Overwrite the ID.qsub files in the subset when mult_jobs==False so if it hits the walltime limit
    in the subset job it returns to the main pipeline when it restarts.
    '''

    main_qsub_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'_%s.qsub'%ID)
    if os.path.exists(main_qsub_file):
        with open(main_qsub_file,'r') as main_qsub_file_obj:
            main_qsub_file_str=main_qsub_file_obj.read()
    else:
        return

    for new_ID,new_oneCalc in list(calc_subset.items()):
        subset_qsub_file= os.path.join(new_oneCalc['_AFLOWPI_FOLDER_'],'_%s.qsub'%new_ID)
        if os.path.exists(subset_qsub_file):
            with open(subset_qsub_file,'w') as subset_qsub_file_obj:
                subset_qsub_file_obj.write(main_qsub_file_str)





def _swap_walltime_logs(main_ID,main_oneCalc,oneCalc,ID):
    '''
    Overwrites the walltime log of the main pipeline job that submitted the subset so the timer is
    correct when the subset calcs return to the main pipeline.
    '''
    subset_walltime_log=AFLOWpi.run._readWalltimeLog(oneCalc,ID)
    AFLOWpi.run._writeWalltimeLog(main_oneCalc,main_ID,subset_walltime_log)
    #move files to local scratch if it's being used for the main pipeline 
    AFLOWpi.prep._to_local_scratch(main_oneCalc,main_ID)    


        
