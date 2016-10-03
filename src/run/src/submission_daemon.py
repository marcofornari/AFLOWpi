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
        qf=re.sub(u'.submit',u'.qsub',submit_file)
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
    for ID,oneCalc in calcs.iteritems():
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

    print 'generating daemon script in %s'%daemon_file_name
    with open(daemon_file_name,'w') as daemon_file_object:
        daemon_file_object.write(daemon_file_string)


    AFLOWpi.run._start_submission_daemon(daemon_file_name)

def _start_submission_daemon(daemon_file_name):
    print 'starting daemon script: %s'%daemon_file_name    
    cur_dir=os.curdir
    os.chdir(os.path.dirname(daemon_file_name))

    subprocess.Popen('nohup python ./submit_daemon.py  &',shell=True)
    os.chdir(cur_dir)
    print 'daemon script started'

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

    for prefix in calc_list_by_chain.keys():
        chain = calc_list_by_chain[prefix]

        for ID in chain.keys():
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
        print 'All calculations in set are finished. Stopping daemon'
        sys.exit(0)
    else:
        return submission_check_list


def _sort_by_chain(calc_list):
    by_chain={}

    for calc_set in calc_list:
        for ID,oneCalc in calc_set.iteritems():        
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
