import AFLOWpi
import os
import json
import logging

class EnvironConfig():
	def __init__(self):
		# read from AFLOW config file and save important things
		self.calcs = None
		self.config = {}

		# the config dictionary here determines the environ specific settings or all calcs. This loads in two different ways.
		# a) calcs dictionary is initialized and has no environ information (aflowpi setup before environ module)
		#    in this case, simply load the calcs, and a default config dictionary
		# b) onecalc contains environ information (run)
		#    in this case, load in the onecalc into self.config

	def init_calcs(self, calcs):
		self.calcs = calcs

	def init_environ(self):
		#workdir = wpre + projectname + '/' + setname + '/'
		self.config['pdict'] = {}
		#self.config['workdir'] = workdir

		# by default set to fa-ionic, maybe change
		self.config['solvent'] = 'water'
		self.config['interface'] = 'ionic'
		self.config['diffuse'] = 'none'

	def from_one_calc(self, oneCalc):
		if "_AFLOW_ENVIRON_" not in oneCalc:
			logging.warning("config.py tried to load oneCalc expecting _AFLOW_ENVIRON_ key, but no key existed")
			return
		self.config = oneCalc["_AFLOW_ENVIRON_"].copy()

	def add_single_loop(self, param, rangelist):
		self.config['pdict']['loopidx'] = 0
		self.config['pdict']['loopvals'] = rangelist
		self.config['loopparam'] = param
		# include a valid python identifier
		self.config['looplabel'] = "loopvals"
		self.config['mode'] = 'loop'

	def edit(self, key, val):
		if 'edit' in self.config:
			self.config['edit'].append([key, val])
		else:
			self.config['edit'] = [[key, val]]

	def set_interface(self, interface):
		self.config['interface'] = interface

	def set_diffuse(self, diffuse):
		self.config['diffuse'] = diffuse

	def update_calcs(self):
		if self.calcs is None:
			logging.error("environ.config.write_to_file called but calcs has not been initialized, this shouldn't happen!")
			return None
		for ID, oneCalc in self.calcs.int_dict.items():
			oneCalc["_AFLOW_ENVIRON_"] = self.config.copy()
			# most importantly save each oneCalc here! (I don't know when they are first created)
			AFLOWpi.prep._saveOneCalc(oneCalc, ID)
		return self.calcs
	
def init_config(calcs):
	environ_config = AFLOWpi.environ.config.EnvironConfig()
	environ_config.init_calcs(calcs)
	environ_config.init_environ()
	return environ_config

def set_params(config):
	astr = ""
	for key, val in list(config.config['pdict'].items()):
		astr += "{} = {}\n".format(str(key), str(val))
	print(('set_params output: {}'.format(astr)))
	return astr

def set_workflow(config, mode, execPrefix, execPostfix):
	scfsingle = '''oneCalc, ID = AFLOWpi.environ._run_environ_single(__submitNodeName__, oneCalc, ID,
		mode="{}", execPrefix="{}", execPostfix="{}")
oneCalc['__execCounter__']+=1
AFLOWpi.prep._saveOneCalc(oneCalc, ID)'''.format(mode, execPrefix, execPostfix)

	astr = ""
	if 'mode' in config.config and config.config['mode'] == 'loop':
		# 1D loop add to script
		param = config.config['loopparam']
		label = config.config['looplabel']
		plist = config.config['pdict']['loopvals']
		for _ in plist:
			astr += ("AFLOWpi.environ.get_environ_input(oneCalc, 'from_config', loopparam='{}', "
					 "loopval={}[loopidx], loopidx=loopidx)\n".format(param, label))
			astr += ("loopidx += 1\n")
			astr += scfsingle + '\n'
			astr += ("shutil.copy(ID+'.out', 'STEP_{:02d}'.format(loopidx))\n")
	else:
		# no loops
		astr += "AFLOWpi.environ.get_environ_input(oneCalc, 'from_config')\n"
		astr += scfsingle + '\n'
	print(('set_workflow output: {}'.format(astr)))
	return astr

