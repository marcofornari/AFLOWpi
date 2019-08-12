import shutil
import json
import os
import logging

class EnvironFile():
	# TODO change defaults to something more general
	# TODO consider changing to a dictionary
	def __init__(self, solvent='vacuum', interface='ionic', diffuse='none'):
		self.edict = {}
		self.set_defaults()
		## EXTERNAL_CHARGES (ignore for now)
		## DIELECTRIC_REGIONS (ignore for now)

		self.set_interface_mode(interface)
		self.set_solvent_mode(solvent)
		self.set_diffuse_mode(diffuse)

	def sanity_check(self, mode):
		pass
		# TODO add sanity checks to make sure input is fine before writing environ.in file

	def set_defaults(self):
		self.edict = {
			## MISC FLAGS
			# diffuse layer check
			'diffuse_layer': False,
			## &ENVIRON namelist
			# environ type is input by default, other options are useful, however, for the purposes of
			# specific tweaking, use input for everything and set special cases via functions
			'environ_type': 'input',
			# verbosity defaults to zero, 1 is useful for the debug file, higher may also be useful
			'verbose': 0,
			# environ_restart defaults to false, if initial guess is good / restart calcs, may want true
			'environ_restart': False,
			# probably useless
			'oldenviron': False,
			# defaults to 0.1, valid for small systems
			'environ_thr': 0.1,
			# not sure what this is for
			'environ_nskip': 1,
			# embedding settings (ignore for now)
			'system_ntyp': 0,
			# embedding dimensionality (ignore for now)
			'system_dim': 0,
			# main axis of embedded system (ignore for now)
			'system_axis': 3,
			# read the electrostatic part? (yes for sscs testing)
			'env_electrostatic': False,
			# static permittivity required if input is set by environ type (default to water)
			'env_static_permittivity': 1.0,
			# for tdddft (ignore for now)
			'env_optical_permittivity': 1.0,
			# surface tension in CGS, required if environ_type set to input
			'env_surface_tension': 0.0,
			# pressure in GPa, required if environ_type set to input
			'env_pressure': 0.0,
			# external charges
			'env_external_charges': 0,
			# dielectric regions
			'env_dielectric_regions': 0,
			# number of electrolyte types
			'env_electrolyte_ntyp': 0,
			# stern correction
			'stern_entropy': 'full',
			# molar concentration, not empty if ionic countercharges are set by env_electrolyte_ntyp
			'cion': [],
			# max molar concentration
			'cionmax': 0.0,
			# mean atomic radius of ionic countercharge
			'rion': 0.0,
			# ionic countercharge, not empty if ionic countercharges are set by env_electrolyte_ntyp
			'zion': [],
			# temperature of electrolyte solution
			'solvent_temperature': 300.0,
			# jellium
			'add_jellium': False,
			# gaussian spread of electrolyte atoms
			'atomicspread': [],
			# linearized electrolyte subroutines
			'electrolyte_linearized': False,
			
			## &BOUNDARY namelist
			# solvent mode sets interface model
			'solvent_mode': 'electronic',
			# if interface function is electronic, below parameters necessary
			'rhomax': 0.005,
			'rhomin': 0.0001,
			'tbeta': 4.8,
			# if interface function is ionic, below parameters necessary
			'alpha': 1.0,
			'softness': 0.5,
			# if interface function is fa-ionic, additional parameters necessary
			'field_awareness': 0.01,
			'charge_asymmetry': 0.0,
			'field_min': 10.0,
			'field_max': 12.0,
			# if interface function is system, below parameters necessary
			'solvent_distance': 1.0,
			'solvent_spread': 0.5,
			# electrolyte
			'electrolyte_mode': 'electronic',
			'electrolyte_rhomax': 5e-3,
			'electrolyte_rhomin': 1e-4,
			'electrolyte_tbeta': 4.8,
			'electrolyte_alpha': 1.0,
			'electrolyte_softness': 0.5,
			'electrolyte_distance': 0.0,
			'electrolyte_spread': 0.5,
			# shape of environment region
			'stype': 1,
			# atomic radii
			'radius_mode': 'uff',
			# manual setting of above
			'solvationrad': [],
			'corespread': [],
			# solvent radius (solvent aware)
			'solvent_radius': 0.0,
			'radial_scale': 2.0,
			'radial_spread': 0.5,
			'filling_threshold': 0.825,
			'filling_spread': 0.02,
			# for interface derivatives
			'boundary_core': 'analytic',
			# see environ input doc
			'ifdtype': 1,
			'nfdpoint': 1,

			## &ELECTROSTATIC namelist
			# electrostatic problem
			'problem': 'poisson',
			# numerical solver
			'solver': 'direct',
			'auxiliary': 'none',
			# accuracy (high = more cycles)
			'tol': 1e-5,
			# mixing parameter
			'mix': 0.5,
			# dimensionality of the simulation cell
			'pbc_dim': 3,
			'pbc_axis': 3,
			'pbc_correction': 'none'
		}

	def set_interface_mode(self, interface):
		# TODO, remove a lot of the assumptions that don't have anything to do with the interface
		if interface == 'electronic':
			self.set_defaults()

			self.edict['solvent_mode'] = 'electronic'
			self.edict['stype'] = 2

			return

		elif interface == 'ionic':
			# set defaults in case they changed (remove in the future?)
			self.set_defaults()

			# set to a standard ionic input (assume vacuum, to be overwritten by set_mode if necessary)
			self.edict['env_electrostatic'] = True

			self.edict['radius_mode'] = 'muff'
			self.edict['solvent_mode'] = 'ionic'
			self.edict['boundary_core'] = 'lowmem'

			self.edict['pbc_correction'] = 'parabolic'
			self.edict['pbc_dim'] = 0
			self.edict['tol'] = 1e-10
			
			return

		elif interface == 'fa-ionic':
			# set defaults in case they changed (remove in the future?)
			self.set_defaults()

			# set to a standard field-aware input (assume water solvent, 
			# to be overwritten by set_mode if necessary)
			self.edict['env_surface_tension'] = 50.0
			self.edict['env_pressure'] = -0.35
			self.edict['env_static_permittivity'] = 78.3
			self.edict['environ_type'] = 'input'
			self.edict['env_electrostatic'] = True

			self.edict['radius_mode'] = 'muff' # seems like uff is not a good idea
			self.edict['solvent_mode'] = 'fa-ionic'
			self.edict['boundary_core'] = 'lowmem'
			self.edict['alpha'] = 1.12
			self.edict['field_awareness'] = 0.12
			self.edict['charge_asymmetry'] = -0.354
			self.edict['field_min'] = 2.0
			self.edict['field_max'] = 3.0

			self.edict['pbc_correction'] = 'parabolic'
			self.edict['pbc_dim'] = 0
			self.edict['tol'] = 1e-10

			return

			# set to a 

	def set_solvent_mode(self, solvent):
		if solvent == 'vacuum':
			# defaults
			self.edict['environ_type'] = 'input'
			self.edict['env_surface_tension'] = 0.0
			self.edict['env_pressure'] = 0.0
			self.edict['env_static_permittivity'] = 1.0

			s = self.edict['solvent_mode']
			if s == 'ionic' or s == 'fa-ionic':
				self.edict['alpha'] = 1.12
			elif s == 'electronic' or s == 'fa-electronic':
				self.edict['rhomin'] = 1e-4
				self.edict['rhomax'] = 5e-3

		elif solvent == 'water':
			self.edict['environ_type'] = 'input'
			self.edict['env_surface_tension'] = 50.0
			self.edict['env_pressure'] = -0.35
			self.edict['env_static_permittivity'] = 78.3

			s = self.edict['solvent_mode']
			if s == 'ionic' or s == 'fa-ionic':
				self.edict['alpha'] = 1.12
			elif s == 'electronic' or s == 'fa-electronic':
				self.edict['rhomin'] = 1e-4
				self.edict['rhomax'] = 5e-3

	def set_diffuse_mode(self, diffuse):
		if diffuse == 'none':
			# defaults
			pass

		elif diffuse == '3d':
			self.edict['diffuse_layer'] = True
			self.edict['env_electrolyte_ntyp'] = 2
			self.edict['zion'] = [1, -1]
			self.edict['cion'] = [1.0, 1.0]
			self.edict['cionmax'] = 10.0
			self.edict['temperature'] = 300.0

			self.edict['electrolyte_mode'] = 'electronic'

			self.edict['pbc_correction'] = 'parabolic'
			self.edict['pbc_dim'] = 0
			self.edict['pbc_axis'] = 3
			self.edict['tol'] = 5e-13
			self.edict['inner_tol'] = 5e-18

			return

	def write_file(self, output, indent=2):
		# ignore most of the options for now, just have workable writefile for SSCS
		s = ''
		for i in range(indent):
			s += ' '
		with open(output, 'w') as f:
			f.write('&ENVIRON\n')
			if self.edict['environ_restart']:
				f.write(s+'environ_restart = .%s.\n'%(str(self.edict['environ_restart']).upper()))
			if self.edict['environ_type'] == 'input':
				# output all input things
				f.write(s+'environ_type = "%s"\n'%(self.edict['environ_type']))
				f.write(s+'env_surface_tension = %f\n'%(self.edict['env_surface_tension']))
				f.write(s+'env_pressure = %f\n'%(self.edict['env_pressure']))
				f.write(s+'env_static_permittivity = %f\n'%(self.edict['env_static_permittivity']))
			if self.edict['env_electrostatic']:
				f.write(s+'env_electrostatic = .%s.\n'%(str(self.edict['env_electrostatic']).upper()))
			if self.edict['environ_thr'] != 0.1:
				f.write(s+'environ_thr = %e\n'%(self.edict['environ_thr']))
			if self.edict['diffuse_layer']:
				f.write(s+'env_electrolyte_ntyp = %d\n'%(self.edict['env_electrolyte_ntyp']))
				for n_electrolyte in range(self.edict['env_electrolyte_ntyp']):
					# TODO more graceful checking.. for now just expect it exists and panic if not
					key = 'zion(%d)'%(n_electrolyte+1)
					f.write(s+'%s = %d\n'%(key, self.edict['zion'][n_electrolyte]))
					key = 'cion(%d)'%(n_electrolyte+1)
					f.write(s+'%s = %f\n'%(key, self.edict['cion'][n_electrolyte]))
				if self.edict['cionmax'] > 0.0:
					f.write(s+'cionmax = %f\n'%(self.edict['cionmax']))
				f.write(s+'temperature = %f\n'%(self.edict['temperature']))
				f.write(s+'electrolyte_linearized = .%s.\n'%(str(self.edict['electrolyte_linearized']).upper()))
			f.write('/\n&BOUNDARY\n')
			if self.edict['radius_mode'] != 'uff':
				f.write(s+'radius_mode = "%s"\n'%(self.edict['radius_mode']))
			if self.edict['solvent_mode'] != 'electronic':
				f.write(s+'solvent_mode = "%s"\n'%(self.edict['solvent_mode']))
			if self.edict['boundary_core'] != 'analytic':
				f.write(s+'boundary_core = "%s"\n'%(self.edict['boundary_core']))
			if self.edict['environ_type'] == 'input':
				if self.edict['solvent_mode'] == 'ionic' or self.edict['solvent_mode'] == 'fa-ionic':
					f.write(s+'alpha = %f\n'%(self.edict['alpha']))
				elif self.edict['solvent_mode'] == 'electronic' or self.edict['solvent_mode'] == 'fa-electronic':
					f.write(s+'rhomin = %e\n'%(self.edict['rhomin']))
					f.write(s+'rhomax = %e\n'%(self.edict['rhomax']))
			if self.edict['stype'] != 1:
				f.write(s+'stype = %d\n'%(self.edict['stype']))
			if self.edict['diffuse_layer']:
				f.write(s+'electrolyte_mode = "%s"\n'%(self.edict['electrolyte_mode']))
			if self.edict['electrolyte_distance'] != 0.0:
				f.write(s+'electrolyte_distance = %f\n'%(self.edict['electrolyte_distance']))
			if self.edict['electrolyte_spread'] != 0.5:
				f.write(s+'electrolyte_spread = %f\n'%(self.edict['electrolyte_spread']))
			if self.edict['electrolyte_rhomin'] != 1e-4:
				f.write(s+'electrolyte_rhomin = %e\n'%(self.edict['electrolyte_rhomin']))
			if self.edict['electrolyte_rhomax'] != 5e-3:
				f.write(s+'electrolyte_rhomax = %e\n'%(self.edict['electrolyte_rhomax']))
			if self.edict['electrolyte_tbeta'] != 4.8:
				f.write(s+'electrolyte_tbeta = %f\n'%(self.edict['electrolyte_tbeta']))
			if self.edict['electrolyte_alpha'] != 1.0:
				f.write(s+'electrolyte_alpha = %f\n'%(self.edict['electrolyte_alpha']))
			if self.edict['solvent_mode'] == 'fa-ionic':
				# add field aware parameters
				f.write(s+'field_awareness = %f\n'%(self.edict['field_awareness']))
				f.write(s+'charge_asymmetry = %f\n'%(self.edict['charge_asymmetry']))
				f.write(s+'field_min = %f\n'%(self.edict['field_min']))
				f.write(s+'field_max = %f\n'%(self.edict['field_max']))
			f.write('/\n')
			if self.edict['env_electrostatic']:
				f.write('&ELECTROSTATIC\n')
				if self.edict['pbc_correction'] != 'none':
					f.write(s+'pbc_correction = "%s"\n'%(self.edict['pbc_correction']))
					f.write(s+'pbc_dim = %d\n'%(self.edict['pbc_dim']))
					f.write(s+'pbc_axis = %d\n'%(self.edict['pbc_axis']))
				if self.edict['tol'] != 1e-5:
					f.write(s+'tol = %e\n'%(self.edict['tol']))
				f.write('/\n')
			elif self.edict['diffuse_layer']:
				f.write('&ELECTROSTATIC\n')
				if self.edict['pbc_correction'] != 'none':
					f.write(s+'pbc_correction = "%s"\n'%(self.edict['pbc_correction']))
					f.write(s+'pbc_dim = %d\n'%(self.edict['pbc_dim']))
					f.write(s+'pbc_axis = %d\n'%(self.edict['pbc_axis']))
				if self.edict['tol'] != 1e-5:
					f.write(s+'tol = %e\n'%(self.edict['tol']))
				f.write(s+'inner_tol = %e\n'%(self.edict['inner_tol']))
				f.write('/\n')

	def edit(self, key, val):
		if not key in self.edict:
			print 'key not found, have you set up the environ file correctly?'
			return
		# TODO sanity check the input
		self.edict[key] = val
		return

def get_environ_input(mode, wdir=None, **kwargs):
	if mode == 'from_file':
		# try and see if there is an existing environ file and if so copy
		# otherwise, load default. Read in environ_config to check if
		# any deviations in settings have been declared
		try:
			shutil.copy(wdir + 'environ.in', os.getcwd() + '/' + 'environ.in')
		except FileNotFoundError:
			get_environ_input('from_config', wdir)
	elif mode == 'from_config':
		# try to read a config file, expect environ.json inside AFLOWpi folder
		aflowdir = os.path.join(os.getcwd(), '../', 'AFLOWpi')
		if not os.path.isdir(aflowdir):
			# unexpected folder structure, send warning
			logging.warning((
				'AFLOWpi folder not found, are you running this function before ' 
				'initializing folder structure?'))
			# for now exit...
			return
		configd = {}
		try:
			with open(aflowdir + '/' + 'environ.json', 'r') as f:
				configd = json.load(f)
		except FileNotFoundError:
			# no config found, just get default
			logging.warning((
				'environ.ini not found, have you initialized the config file?'))
			# for now, just load a default
			get_environ_input('default', wdir)

		# read dictionary and parse contents into a template environ file
		interface = configd['interface']
		solvent = configd['solvent']
		diffuse = configd['diffuse']
		efile = EnvironFile(solvent=solvent, interface=interface, diffuse=diffuse)
		# edit depending on config file
		if 'edit' in configd and configd['edit']:
			for edit in configd['edit']:
				efile.edit(edit[0], edit[1])

		if not kwargs:
			# no loop, just write the environ file
			efile.write_file(os.getcwd() + '/' + 'environ.in')
			return

		print kwargs
		if 'loopval' in kwargs and 'loopparam' in kwargs:
			# part of a loop, thus do a substitution based on the param fed in
			print 'updating'
			efile.edict[kwargs['loopparam']] = kwargs['loopval']
			efile.write_file(os.getcwd() + '/' + 'ENVIRON_%02d'%(kwargs['loopidx']+1))
		efile.write_file(os.getcwd() + '/' + 'environ.in')

if __name__ == '__main__':
	e = EnvironFile(interface='electronic', solvent='water', diffuse='3d')
	e.write_file('environ.in')
