from .loop import loop
from .thermometry import thermometry
from ..metadata import metadata


import urllib2
import scipy.io #for loading .mat file
import os
import numpy as np



import tables
import matplotlib as mpl
mpl.use("pgf")
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": [
         r"\usepackage[utf8x]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         #r"\usepackage{cmbright}",
         ]
}

# pgf_with_latex = {                      # setup matplotlib to use latex for output
#     "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
#     "text.usetex": True,                # use LaTeX to write all text
#     "font.family": "serif",
#     "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
#     "font.sans-serif": [],
#     "font.monospace": [],
#     "axes.labelsize": 10,               # LaTeX default is 10pt font.
#     "text.fontsize": 10,
#     "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
#     "xtick.labelsize": 8,
#     "ytick.labelsize": 8,
#     #"figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
#     "pgf.preamble": [
#         r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
#         r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
#         ]
#     }


mpl.rcParams.update(pgf_with_pdflatex)
import matplotlib.pyplot as plt
# plt.switch_backend('pgf')



import matplotlib.gridspec as gridspec




import datetime
from scipy.optimize import minimize, curve_fit, leastsq# , root,newton_krylov, anderson
from scipy.interpolate import interp1d
from scipy import constants

import numpy.ma as ma
import sys # for status percentage

import platform
mysys = platform.system()

import warnings #trying to get a warning every time rather than just the first time.
warnings.filterwarnings('always')

database_location = 'Data' + os.sep + 'My_Data_Library.h5'
Working_Dir = os.getcwd()
Plots_Dir = '/Users/miguel_daal/Documents/Projects/Thesis/Thesis/chap6/images/plots'
if os.path.exists(Plots_Dir) == False:
	print('Speficied plots directory does not exist... Using current directory')
	Plots_Dir = Working_Dir


try:
	execfile('KIPs_Access.txt')
except:
	print('KIPs_Access.txt not found. Create this file to download data.')
	# remote access file must have these two lines
	# username = _________  # e.g. 'johndoe'
	# password = _________  # e.g. '294hr5'
	##############
