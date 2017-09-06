import matplotlib.pyplot as plt

def plot_loop(metadata, loop, aspect='auto', show = True):
	''' Plots currently selected complex transmission in the I,Q plane. Reutrns a tuple, (fig, ax, line),
	where  fig is the figure object, ax is the axes object and line is the line object for the plotted data.

	aspect='equal' makes circles round, aspect='auto' fills the figure space.

	*Must have a loop picked in order to use this function.*
	'''
	try: 
		z = loop.z
	except:
		print("Data not available. Maybe you forgot to load it.")
		return
	fig = plt.figure( figsize=(6.5, 6.5), dpi=100)
	ax = fig.add_subplot(111,aspect=aspect)
	line, = ax.plot(z.real,z.imag,'bo')
	ax.set_xlabel(r'$\Re[S_{21}(f)]$')
	ax.set_ylabel(r'$\Im[S_{21}(f)]$')
	ax.yaxis.labelpad = -2
	ax.set_title('Run: {0}; Sensor: {1}; Ground: {2}; Record Date: {3}'.format(metadata.Run, metadata.Sensor, metadata.Ground_Plane, metadata.Time_Created),fontsize=10)
	

	if show == True:
		plt.show()
	return  (fig, ax, line)
