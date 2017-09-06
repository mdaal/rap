
def pick_loop(loop,Sweep_Array,index):
	'''Use this function to pick the current loop/transmission data from withing the Sweep_Array. 
	Index is the indes number of sweep/loop to be slected as the current loop.'''
	loop.index = index 
	#self.loop.normalization = None
	loop.z = Sweep_Array[index]['S21']
	loop.freq = Sweep_Array[index]['Frequencies']