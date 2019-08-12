class Environment():
	"""
	ENVIRONMENT, container for custom environments to be read in

	Currently has no functionality but exists to be extended
	"""
	# TODO add persistence so that environments can be stored, perhaps
	# in a database
	def __init__(self, name):
		self.name = name
		# flags that show whether parameterization is set
		self.pmted_system = False
		self.pmted_ionic = False
		self.pmted_electronic = False
		# parameterization variables
		self.pressure = 0.0
		self.tension = 0.0
		self.alpha = 1.12
		self.rhomax = 5e-3
		self.rhomin = 1e-4

	def set_parameterization(self, pressure=None, tension=None,
			alpha=None, rhomax=None, rhomin=None):
		system = False
		ionic = False
		electronic = False
		# if None just do nothing
		if pressure is not None:
			self.pressure = pressure
			system = True
		if tension is not None:
			self.tension = tension
			system = True
		if alpha is not None:
			self.alpha = alpha
			ionic = True
		if rhomax is not None:
			self.rhomax = rhomax
			electronic = True
		if rhomin is not None:
			self.rhomin = rhomin
			electronic = True
		if electronic:
			self.pmted_electronic = True
		if ionic:
			self.pmted_ionic = True
		if system:
			self.pmted_system = True

