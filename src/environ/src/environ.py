import AFLOWpi
import os
import logging

def _execheck():
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
	return pwx_exec_loc, symlink

def setup_relax(calcs, config=None, environmode='from_file'):
	"""SETUP_RELAX
	this function should only need to worry about making sure
	calculation is a relax and that the environ.in file is
	correctly copied over
	"""

	print "entering SETUP"
	if config is not None:
		print "config set to not None"

	# TODO remove commented code once testing is complete

	#step_counter=0

	#AFLOWpi.prep.addToAll_(calcs,'PREPROCESSING',"""oneCalc,ID=AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"scf"')""")

	pwx, symlink = _execheck()

	AFLOWpi.prep.totree(pwx, calcs, rename=None, symlink=symlink)

	#step_counter+=1

	#command='''if oneCalc["__execCounter__"]<=%s:
	#oneCalc,ID=AFLOWpi.environ._run_environ_relax(__submitNodeName__,oneCalc,ID)
	#oneCalc['__execCounter__']+=1
	#AFLOWpi.prep._saveOneCalc(oneCalc,ID)'''%(step_counter)
	
	execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
	execPostfix=AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
	runcommand='''oneCalc, ID = AFLOWpi.environ._run_environ_relax(__submitNodeName__, oneCalc, ID, execPrefix="%s",
		execPostfix="%s")
oneCalc['__execCounter__']+=1
AFLOWpi.prep._saveOneCalc(oneCalc, ID)'''%(execPrefix, execPostfix)

	working_directory = os.getcwd() + "/"
	param_pre = """wdir = '%s'"""%(working_directory)
	if config is not None:
		param_pre = AFLOWpi.environ.set_params(config)
	if config is None:
		environ_pre = """AFLOWpi.environ.get_environ_input('%s', wdir)"""%(environmode)
	else:
		environ_pre = """import shutil"""

	AFLOWpi.prep.addToAll_(calcs, 'PREPROCESSING', param_pre)
	AFLOWpi.prep.addToAll_(calcs, 'PREPROCESSING', environ_pre)

	if config is not None:
		runcommand = AFLOWpi.environ.set_workflow(config, 'relax', execPrefix, execPostfix)
	AFLOWpi.prep.addToAll_(calcs, 'RUN', runcommand)

def setup_scf(calcs, config=None, environmode='from_file'):
	"""SETUP_SCF
	this function needs to make sure the calculation is scf
	and handle a change of environ file before each step of the
	calculation is made.

	Since this is technically not a multi-step calculation, maybe just
	have everything run iteratively from one python file... (IE do everything from
	RUN phase, test for this.

	TODO: merge with relax (no reason to have a setup for each)
	"""

	print "entering SETUP"
	if config is not None:
		print "config set to not None"

	pwx, symlink = _execheck()

	AFLOWpi.prep.totree(pwx, calcs, rename=None, symlink=symlink)

	execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
	execPostfix=AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
	runcommand='''oneCalc, ID = AFLOWpi.environ._run_environ_scf(__submitNodeName__, oneCalc, ID, execPrefix="%s",
		execPostfix="%s")
oneCalc['__execCounter__']+=1
AFLOWpi.prep._saveOneCalc(oneCalc, ID)'''%(execPrefix, execPostfix)

	working_directory = os.getcwd() + "/"
	param_pre = """wdir = '%s'"""%(working_directory)
	if config is not None:
		param_pre = AFLOWpi.environ.set_params(config)
	if config is None:
		environ_pre = """AFLOWpi.environ.get_environ_input('%s', wdir)"""%(environmode)
	else:
		environ_pre = """import shutil"""

	AFLOWpi.prep.addToAll_(calcs, 'PREPROCESSING', param_pre)
	AFLOWpi.prep.addToAll_(calcs, 'PREPROCESSING', environ_pre)

	if config is not None:
		runcommand = AFLOWpi.environ.set_workflow(config, 'scf', execPrefix, execPostfix)
	AFLOWpi.prep.addToAll_(calcs, 'RUN', runcommand)

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
		execPostfix = ""
		execPostfix += " --environ "

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
		AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix, 
                execPostfix=execPostfix, engine='espresso', calcType='scf', executable=None, exit_on_error=False)

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

		AFLOWpi.run._oneRun(__submitNodeName__, environ_scf_calc, environ_scf_ID, execPrefix=execPrefix,
			execPostfix=execPostfix, engine='espresso', calcType='scf', executable=None, exit_on_error=False)

		oneCalc['__runList__'].append('environ_scf')
		AFLOWpi.prep._saveOneCalc(oneCalc,ID)	

		return oneCalc,ID

##################################################################################################################

def _run_environ_relax(__submitNodeName__, oneCalc, ID, execPrefix, execPostfix):
	oneCalcID = ID
	engine = ''
	config = None

	if "__runList__" not in oneCalc.keys():
		oneCalc["__runList__"] = []
		config = None

	if config != None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug("forced config %s" % config)
	else:
		try:
			config = AFLOWpi.prep._getConfigFile()
			AFLOWpi.prep._forceGlobalConfigFile(config)
		except Exception, e:
			AFLOWpi.run._fancy_error_log(e)

	if AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix") != "":
		execPostfix = AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix")
	else:
		execPostfix = ""
		execPostfix += " --environ"

	if AFLOWpi.prep._ConfigSectionMap("run", "engine") == "":
		engine = AFLOWpi.prep._ConfigSectionMap("run", "engine")
	else:
		engine = "espresso"

	oneCalc["_AFLOWPI_CONFIG_"] = config

	# RELAX
	oneCalc, ID = AFLOWpi.environ._setup_environ_relax(oneCalc, ID)

	if "environ_relax" not in oneCalc["__runList__"]:
		AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix,
				execPostfix=execPostfix, engine=engine, calcType="scf",
				executable=None, exit_on_error=False)

		oneCalc["__runList__"].append("environ_relax")
		AFLOWpi.prep._saveOneCalc(oneCalc, ID)

		return oneCalc, ID
	else:
		# ignore above check for now because I don't understand why it is defaulting
		AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix,
				execPostfix=execPostfix, engine=engine, calcType="scf",
				executable=None, exit_on_error=False)

		oneCalc["__runList__"].append("environ_scf")
		AFLOWpi.prep._saveOneCalc(oneCalc, ID)

		return oneCalc, ID

def _run_environ_scf(__submitNodeName__, oneCalc, ID, execPrefix, execPostfix):
	oneCalcID = ID
	engine = ''
	config = None

	if "__runList__" not in oneCalc.keys():
		oneCalc["__runList__"] = []
		config = None

	if config != None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug("forced config %s" % config)
	else:
		try:
			config = AFLOWpi.prep._getConfigFile()
			AFLOWpi.prep._forceGlobalConfigFile(config)
		except Exception, e:
			AFLOWpi.run._fancy_error_log(e)

	if AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix") != "":
		execPostfix = AFLOWpi.prep._ConfigSectionMap("run", "exec_postfix")
	else:
		execPostfix = ""
		execPostfix += " --environ"

	if AFLOWpi.prep._ConfigSectionMap("run", "engine") == "":
		engine = AFLOWpi.prep._ConfigSectionMap("run", "engine")
	else:
		engine = "espresso"

	oneCalc["_AFLOWPI_CONFIG_"] = config

	# RELAX
	oneCalc, ID = AFLOWpi.environ._setup_environ_scf(oneCalc, ID)

	if "environ_scf" not in oneCalc["__runList__"]:
		AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix,
				execPostfix=execPostfix, engine=engine, calcType="scf",
				executable=None, exit_on_error=False)

		oneCalc["__runList__"].append("environ_scf")
		AFLOWpi.prep._saveOneCalc(oneCalc, ID)

		return oneCalc, ID
	else:
		# ignore above check for now because I don't understand why it is defaulting
		AFLOWpi.run._oneRun(__submitNodeName__, oneCalc, ID, execPrefix=execPrefix,
				execPostfix=execPostfix, engine=engine, calcType="scf",
				executable=None, exit_on_error=False)

		oneCalc["__runList__"].append("environ_scf")
		AFLOWpi.prep._saveOneCalc(oneCalc, ID)

		return oneCalc, ID

def _setup_environ_scf(oneCalc,ID):

	onecalc, ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"scf"')

	return oneCalc, ID

def _setup_environ_relax2scf(oneCalc,ID):

	environ_scf_ID='%s_environ_scf' % ID
	environ_scf_oneCalc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],ID)                    

	# CHANGE THE RELAX TO SCF AND SAVE FILE TO DISK
	environ_scf_oneCalc,environ_scf_ID=AFLOWpi.prep._modifyNamelistPW(environ_scf_oneCalc,environ_scf_ID,'&control','calculation','"scf"')    

	return environ_scf_oneCalc,environ_scf_ID


def _setup_environ_relax(oneCalc,ID):

	# CHANGE TO RELAX AND SAVE FILE TO DISK
	oneCalc,ID=AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','calculation','"relax"')    

	return oneCalc,ID


