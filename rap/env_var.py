import platform
import os

class env_var():
	def __init__(self, **kword):
		self.mysys = platform.system()

		if  'database_location' in kword.keys():
			self.database_location = kword['database_location']
		else:
			self.database_location = os.path.expanduser('~') + os.sep + 'Resonator_Data' + os.sep + 'My_Data_Library.h5'
