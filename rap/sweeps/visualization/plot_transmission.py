import matplotlib.pyplot as plt
import numpy as np

def plot_transmission(metadata, loop,  show = True):
	''' Plots currently selected complex transmission in dB as a function of frequency. Reutrns a tuple, (fig, ax, line),
	where fig is the figure object, ax is the axes object and line is the line object for the plotted data.

	*Must have a loop picked in order to use this function.*
	'''
	try: 
		z = loop.z
		freq = loop.freq
	except:
		print("Data not available. You probably forgot to load it.")
		return
	plt.rcParams["axes.titlesize"] = 10
	fig = plt.figure( figsize=(8, 6), dpi=100)
	ax = fig.add_subplot(111)
	line = ax.plot(freq,20*np.log10(abs(z)),'b-',)
	ax.set_xlabel('Frequency [Hz]')
	ax.set_ylabel('$20*Log_{10}[|S_{21}|]$ [dB]')

	ax.set_title('Run: {0}; Sensor: {1}; Ground: {2}; Record Date: {3}'.format(metadata.Run, metadata.Sensor, metadata.Ground_Plane, metadata.Time_Created))
	if show == True:
		plt.show()
	return  (fig, ax, line)
