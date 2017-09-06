import numpy as np

def normalize_loop(loop, base, offset):
	''' normalize loop so that mag(S21)< 1. determine normalization by averaging np.abs(S21[base:offset]).mean()
	return normalization'''
	S21 = loop.z
	f = loop.freq	
	
	normalization = np.abs(S21[base:offset]).mean() # consider using medium()?
	loop.normalization = normalization
	S21_normalized = S21/normalization
	loop.z = S21_normalized

	return normalization 