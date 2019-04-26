import AFLOWpi
import os

def _setup_environ(calcs):

        step_counter=0

        AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""oneCalc,ID=AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"scf"')""")


        pwx_dir=AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
            symlink=True
        else:
            symlink=False
        pwx_exec_loc = os.path.join(pwx_dir,'pw.x')
        if not os.path.exists(pwx_exec_loc):
            logging.error('ERROR: engine executables not found in %s please check your config file. EXITING' % pwx_dir)
            print 'ERROR: engine executables not found in %s please check your config file EXITING' % pwx_dir
            raise SystemExit

        AFLOWpi.prep.totree(pwx_exec_loc, calcs,rename=None,symlink=symlink)



        command='''if oneCalc["__execCounter__"]<=%s:
     oneCalc,ID=AFLOWpi.environ._run_environ_iterative(__submitNodeName__,oneCalc,ID)
     oneCalc['__execCounter__']+=1
     AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''%(step_counter)

        AFLOWpi.prep.addToAll_(calcs,'RUN',command) 

        step_counter+=1



def _run_environ_iterative(__submitNodeName__,oneCalc,ID):

##################################################################################################################
# probably don't have to worry about this stuff
##################################################################################################################
	execPrefix = ''
	execPostfix = ''
        oneCalcID = ID

        if '__runList__' not in oneCalc.keys():
            oneCalc['__runList__']=[]

        config=None
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

        execPostfix+=" -environ "


	if AFLOWpi.prep._ConfigSectionMap('run','engine') == '':
		engine = AFLOWpi.prep._ConfigSectionMap('run','engine')
	else:
		engine = 'espresso'

	oneCalc['_AFLOWPI_CONFIG_']=config


##################################################################################################################
# do relax
##################################################################################################################

        oneCalc,ID = AFLOWpi.environ._setup_environ_relax(oneCalc,ID)

        if 'environ_relax' not in oneCalc['__runList__']:

            AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)

            oneCalc['__runList__'].append('environ_relax')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            
            environ_scf_calc,environ_scf_ID= AFLOWpi.environ._setup_environ_scf(oneCalc,ID)


        else:
            '''if we are restarting from a job killed from going walltime 
            try to load environ_scf_ID and if we can't then just make a new one'''
            try:
                environ_scf_ID='%s_environ_scf' % ID
                environ_scf_calc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],environ_scf_ID)                
                '''we have to make sure nscf step has the correct walltime and start time if it's a restart'''
                environ_scf_calc['__walltime_dict__']=oneCalc['__walltime_dict__']
            except Exception,e:
                try:
                    # setup the scf after the relax
                    environ_scf_calc,environ_scf_ID= AFLOWpi.environ._setup_environ_scf(oneCalc,ID)
                                                                       	

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)


##################################################################################################################
# do scf
##################################################################################################################

        if 'environ_scf' not in oneCalc['__runList__']:


            AFLOWpi.run._oneRun(__submitNodeName__,environ_scf_calc,environ_scf_ID,execPrefix=execPrefix,
                                execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)



            oneCalc['__runList__'].append('environ_scf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)	


        return oneCalc,ID

##################################################################################################################



def _setup_environ_scf(oneCalc,ID):

    environ_scf_ID='%s_environ_scf' % ID
    environ_scf_oneCalc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],ID)                    

    # CHANGE THE RELAX TO SCF AND SAVE FILE TO DISK
    environ_scf_oneCalc,environ_scf_ID=AFLOWpi.prep._modifyNamelistPW(environ_scf_oneCalc,environ_scf_ID,'&control','calculation','"scf"')    

    return environ_scf_oneCalc,environ_scf_ID


def _setup_environ_relax(oneCalc,ID):

    # CHANGE TO RELAX AND SAVE FILE TO DISK
    oneCalc,ID=AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"relax"')    

    return oneCalc,ID


