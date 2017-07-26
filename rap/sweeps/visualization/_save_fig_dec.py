def _save_fig_dec(self, fig, name, Use_Date = False, Make_PGF = True):
		os.chdir(Plots_Dir)
		if self.metadata.Run is not None:
			name = self.metadata.Run+ '_'+ name  
		if Use_Date:
			name = name + '_'+ datetime.date.today().strftime("%Y%m%d")

		fig.savefig(name, dpi=300, transparency  = True, bbox_inches='tight')#Title.replace('\n','_').replace(' ','_')+date
		if Make_PGF:
			#cur_backend = mpl.get_backend()
			#plt.switch_backend('pgf')
			name = name + '.pgf'
			plt.savefig(name, bbox_inches = 'tight', transparancy = True) #
			#plt.switch_backend(cur_backend)
		os.chdir(Working_Dir)

