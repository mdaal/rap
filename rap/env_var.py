import platform
import os
import matplotlib as mpl

class env_var():
	def __init__(self, **kword):
		self.mysys = platform.system()

		if  'database_location' in kword.keys():
			self.database_location = kword['database_location']
		else:
			self.database_location = os.path.expanduser('~') + os.sep + 'Resonator_Data' + os.sep + 'My_Data_Library.h5'



		pgf_with_pdflatex = {
			"interactive": True,
			"pgf.texsystem": "pdflatex",
			"pgf.preamble": [
				r"\usepackage[utf8x]{inputenc}",
				r"\usepackage[T1]{fontenc}",
				#r"\usepackage{cmbright}",
				]
			}

		mpl.rcParams.update(pgf_with_pdflatex)
		#self.imported_data_path = None not necessary: using metadate.Data_Source = None
