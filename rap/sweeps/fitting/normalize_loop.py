def normalize_loop(self, base = 0, offset = 5):
	''' normalize loop so that mag(S21)< 1. determine normalization by averaging np.abs(S21[base:offset]).mean()
	return normalization'''
	S21 = self.loop.z
	f= self.loop.freq	
	
	normalization = np.abs(S21[base:offset]).mean() # consider using medium()?
	self.loop.normalization = normalization
	S21_normalized = S21/normalization
	self.loop.z = S21_normalized

	return normalization 