def downsample_loop(N):
	''' Reduce number of loop/freq data point by every Nth point and discarding all others'''
	loop.z = loop.z[0:-1:N]
	loop.freq =loop.freq[0:-1:N]
	loop_downsample=[loop.z,loop.freq]
	return loop_downsample
