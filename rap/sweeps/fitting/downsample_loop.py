def downsample_loop(self,N):
	''' Reduce number of loop/freq data point by every Nth point and discarding all others'''
	self.loop.z = self.loop.z[0:-1:N]
	self.loop.freq = self.loop.freq[0:-1:N]
