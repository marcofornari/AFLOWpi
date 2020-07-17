import AFLOWpi
import os
import logging

def _execheck():
	pwx_dir=AFLOWpi.prep._ConfigSectionMap('prep', 'engine_dir')
	if AFLOWpi.prep._ConfigSectionMap('prep', 'copy_execs').lower() == 'false':
		symlink = True
	else:
		symlink = False
	pwx_exec_loc = os.path.join(pwx_dir, 'pw.x')
	if not os.path.exists(pwx_exec_loc):
		logging.error('ERROR: engine executables not found in {} please check your config file. EXITING'.format(pwx_dir))
		print(('ERROR: engine executables not found in {} please check your config file EXITING'.format(pwx_dir)))
		raise SystemExit
	return pwx_exec_loc, symlink

def setup_environ(calcs, workflow, config=None):
	"""MAIN setup function

	this function is called by the user indirectly in order to write the submission scripts with some
	predetermined workflow. Currently support the basic options, relax and scf.

	Args:
		workflow (str): the flag that determines a desired standard workflow
		config (class:`AFLOWpi.environ.config.EnvironConfig`): an optional config file 
			with non-standard workflow information

	Returns:
		None
	"""
	environmode = "from_file"
	if config is None:
		environmode = "from_config"

	print("entering SETUP, config={}".format(config))
	if config is not None:
		print("config set to not None")
		config.init_calcs(calcs)

	# TODO remove commented code once testing is complete

	#step_counter=0

	#AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""oneCalc,ID=AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"scf"')""")

	pwx, symlink = _execheck()

	AFLOWpi.prep.totree(pwx, calcs, rename=None, symlink=symlink)

	#step_counter+=1

	# command='''if oneCalc["__execCounter__"]<=%s:
	# oneCalc,ID=AFLOWpi.environ._run_environ_relax(__submitNodeName__,oneCalc,ID)
	# oneCalc['__execCounter__']+=1
	# AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''%(step_counter)
	
	runcommand = '''oneCalc, ID = AFLOWpi.environ._run_environ_single(__submitNodeName__, oneCalc, ID, "{}")
oneCalc['__execCounter__']+=1
AFLOWpi.prep._saveOneCalc(oneCalc, ID)'''.format(workflow)

	working_directory = os.getcwd() + "/"
	param_pre = """wdir = '{}'""".format(working_directory)
	if config is not None:
		param_pre = AFLOWpi.environ.set_params(config)
	if config is None:
		# oneCalc should be initialized by the script already
		environ_pre = """AFLOWpi.environ.get_environ_input(oneCalc, '{}', wdir)""".format(environmode)
	else:
		environ_pre = """import shutil"""

	AFLOWpi.prep.addToAll_(calcs, 'PREPROCESSING', param_pre)
	AFLOWpi.prep.addToAll_(calcs, 'PREPROCESSING', environ_pre)

	if config is not None:
		runcommand = AFLOWpi.environ.set_workflow(config, workflow)
	AFLOWpi.prep.addToAll_(calcs, 'RUN', runcommand)

def _run_environ_iterative(__submitNodeName__, oneCalc, ID):
	"""INTERNAL iterative function

	This function is untested, and reserved since it contains potentially useful info on how
	to construct particular workflows. Intended to run a relax and scf sequentially

	Args:
		__submitNodeName__ ([type]): [description]
		oneCalc (dict): information specific to the submission
		ID (str): unique job identifier

	Returns:
		tuple(dict, str): returns updated versions of input parameters
	"""

	# probably don't have to worry about this stuff
	execPrefix = ''
	execPostfix = ''
	oneCalcID = ID

	if '__runList__' not in list(oneCalc.keys()):
		oneCalc['__runList__'] = []
		config=None

	if config is not None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug('forced config {}'.format(config))
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

	if AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix") != '':
		execPostfix = AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix")
	else:
		execPostfix = ""
		execPostfix += " --environ "

	if AFLOWpi.prep._ConfigSectionMap('run', 'engine') == '':
		engine = AFLOWpi.prep._ConfigSectionMap('run', 'engine')
	else:
		engine = 'espresso'

	oneCalc['_AFLOWPI_CONFIG_'] = config

	oneCalc, ID = AFLOWpi.environ._setup_environ_relax(oneCalc, ID)

	if 'environ_relax' not in oneCalc['__runList__']:
		AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix, 
                execPostfix=execPostfix, engine='espresso', calcType='scf', executable=None, exit_on_error=False)

		oneCalc['__runList__'].append('environ_relax')
		AFLOWpi.prep._saveOneCalc(oneCalc,ID)
        
		environ_scf_calc, environ_scf_ID = AFLOWpi.environ._setup_environ_scf(oneCalc, ID)

	else:
		'''if we are restarting from a job killed from going walltime 
		   try to load environ_scf_ID and if we can't then just make a new one'''
		try:
			environ_scf_ID='{}_environ_scf'.format(ID)
			environ_scf_calc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],environ_scf_ID)                
			'''we have to make sure nscf step has the correct walltime and start time if it's a restart'''
			environ_scf_calc['__walltime_dict__']=oneCalc['__walltime_dict__']
		except Exception as e:
			try:
			    # setup the scf after the relax
				environ_scf_calc,environ_scf_ID= AFLOWpi.environ._setup_environ_scf(oneCalc,ID)
                                                                
			except Exception as e:
				AFLOWpi.run._fancy_error_log(e)

	if 'environ_scf' not in oneCalc['__runList__']:

		AFLOWpi.run._oneRun(__submitNodeName__, environ_scf_calc, environ_scf_ID, execPrefix=execPrefix,
			execPostfix=execPostfix, engine='espresso', calcType='scf', executable=None, exit_on_error=False)

		oneCalc['__runList__'].append('environ_scf')
		AFLOWpi.prep._saveOneCalc(oneCalc, ID)	

	return oneCalc, ID

def _run_environ_single(__submitNodeName__, oneCalc, ID, mode):
	"""INTERNAL single function

	This function runs a single environ calculation, determined by the mode param, which
	currently switches between scf and relax

	Args:
		__submitNodeName__ ([type]): [description]
		oneCalc (dict): information specific to the submission
		ID (str): unique job identifier
		mode (str): job type

	Returns:
		tuple(dict, str): returns updated versions of input parameters
	"""
	oneCalcID = ID
	engine = ''
	config = None

	if "__runList__" not in list(oneCalc.keys()):
		oneCalc["__runList__"] = []
		config = None

	if config is not None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug("forced config {}".format(config))
	else:
		try:
			config = AFLOWpi.prep._getConfigFile()
			AFLOWpi.prep._forceGlobalConfigFile(config)
		except Exception as e:
			AFLOWpi.run._fancy_error_log(e)

	if AFLOWpi.prep._ConfigSectionMap("run", "exec_prefix") != '':
		execPrefix=AFLOWpi.prep._ConfigSectionMap("run", "exec_prefix")
	else:
		execPrefix=''

	if AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix") != "":
		execPostfix = AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix")
	else:
		execPostfix = " --environ"

	if AFLOWpi.prep._ConfigSectionMap("run", "engine") == "":
		engine = AFLOWpi.prep._ConfigSectionMap("run", "engine")
	else:
		engine = "espresso"

	oneCalc["_AFLOWPI_CONFIG_"] = config

	if mode == "environ_scf":
		oneCalc, ID = AFLOWpi.environ._setup_environ_scf(oneCalc, ID)
	elif mode == "environ_relax":
		oneCalc, ID = AFLOWpi.environ._setup_environ_relax(oneCalc, ID)

	# TODO check what this does exactly
	AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix,
			execPostfix=execPostfix, engine=engine, calcType="scf",
			executable=None, exit_on_error=False)

	oneCalc["__runList__"].append(mode)
	AFLOWpi.prep._saveOneCalc(oneCalc, ID)

	return oneCalc, ID

def _setup_environ_scf(oneCalc, ID):
	"""INTERNAL setup function

	Not meant to be called directly, ensures the correct parameter entry in the pw input
	file for an scf calculation

	Args:
		oneCalc (dict): information specific to the submission
		ID (str): unique job identifier

	Returns:
		tuple(dict, str): returns updated versions of input parameters
	"""
	oneCalc, ID = AFLOWpi.prep._modifyNamelistPW(oneCalc, ID, '&control', 'calculation', '"scf"')
	return oneCalc, ID

def _setup_environ_relax2scf(oneCalc, ID):
	"""INTERNAL setup function

	Not meant to be called directly, ensures the correct parameter entry in the pw input
	file for an scf calculation

	Args:
		oneCalc (dict): information specific to the submission
		ID (str): unique job identifier

	Returns:
		tuple(dict, str): returns updated versions of input parameters
	"""

	environ_scf_ID = "{}_environ_scf".format(ID)
	environ_scf_oneCalc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'], ID)                    

	# CHANGE THE RELAX TO SCF AND SAVE FILE TO DISK
	environ_scf_oneCalc, environ_scf_ID = AFLOWpi.prep._modifyNamelistPW(
		environ_scf_oneCalc, environ_scf_ID, '&control', 'calculation', '"scf"')    

	return environ_scf_oneCalc, environ_scf_ID


def _setup_environ_relax(oneCalc,ID):
	"""INTERNAL setup function

	Not meant to be called directly, ensures the correct parameter entry in the pw input
	file for a relax calculation

	Args:
		oneCalc (dict): information specific to the submission
		ID (str): unique job identifier

	Returns:
		tuple(dict, str): returns updated versions of input parameters
	"""
	oneCalc, ID = AFLOWpi.prep._modifyNamelistPW(oneCalc, ID, '&control', 'calculation', '"relax"')    
	return oneCalc, ID


