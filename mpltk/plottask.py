
	
	def __plot_light_curve(self,lc_info,t0 = 2457000,xylim=False):
	
		def jd2iso(jd):
			from astropy.time import Time
			t = Time(jd,scale='utc',format='jd')
			t_iso = t.iso
			t_iso_short = t_iso.split()[0]
		
			return t_iso_short
	
	
	        colors = ['g','b','r','y','c','k','m']
		symbols = ['o','s','D','*','+','p','H']
	#        ax = plt.subplot(111)

		
		ax = host_subplot(111, axes_class=AA.Axes)
	
		target_name = lc_info['sn']
		lc_data_dir = lc_info['lc_spec_dir']
		lc_filenames = lc_info['lc']
	
		cm = {}
		ci = 0 
		mi = 0	
	

		t_min = None
		t_max = None
		m_min = None
		m_max = None
		


	        for i,lc_file in enumerate(lc_filenames):
	
		    
			instrument_band = lc_file.split('_')[1].split('.')[0]
	        	instrument  = instrument_band.split('-')[0]
			flt = instrument_band.split('-')[1]
			if flt == 'uu':
				flt = 'U'
			if flt == 'bb':
	                	flt = 'B'
			if flt == 'vv':
	                	flt = 'V'
	
			print flt
	
			if instrument not in cm.keys():
				cm[instrument] = mi
				mi = mi+1
			if flt not in cm.keys():
				
				cm[flt] = ci
				ci = ci+1
	 
		
	
			lc = np.loadtxt(os.path.join(lc_data_dir,lc_file))
			
			try:
				t     = lc[:,0]-t0
				mag   = lc[:,1]
				magerr= lc[:,2]
			except:
				t = lc[0] - t0
				mag = lc[1]
				magerr = lc[2]
			
			t_min_this = np.min(t)
			t_max_this = np.max(t)
			m_min_this = np.min(mag)
			m_max_this = np.max(mag)

			if t_min is None or t_min_this < t_min:
				t_min = t_min_this
			if t_max is None or t_max_this > t_max:
				t_max = t_max_this             			
			if m_min is None or m_min_this < m_min:
				m_min = m_min_this
			if m_max is None or m_max_this > m_max:
				m_max = m_max_this

			plt_hand = self.__plot_light_curve_single(ax,t,mag,err = magerr,color = colors[np.mod(cm[flt],7)],marker = symbols[np.mod(cm[instrument],7)])
	
	
		print t_min,t_max,m_min,m_max
		
		if xylim:	
			plt.xlim([t_min-0.8,t_max+0.8])
			plt.ylim([m_min-0.5,m_max+0.5])
	
		xticks = plt.xticks()
		ax2 = ax.twin() # ax2 is responsible for "top" axis and "right" axis
		ax2.set_xticks(xticks[0])
	
		ax2ticklabels = [jd2iso(xticklabel +t0 ) for xticklabel in xticks[0]]
		
		ax2ticklabels_new = []
		suffix_prior = ''
		for ticklabel in ax2ticklabels:
		 	ticklabel_coms = ticklabel.split('-')
			suffix =ticklabel_coms[0]
			if suffix == suffix_prior:
				ax2ticklabels_new.append(ticklabel_coms[1]+'-'+ticklabel_coms[2])
			else:
				suffix_prior = ticklabel_coms[0]
				ax2ticklabels_new.append(ticklabel_coms[0]+'-'+ticklabel_coms[1]+'-'+ticklabel_coms[2])
		
		ax2.set_xticklabels(ax2ticklabels_new)
		
		ax2.axis["right"].major_ticklabels.set_visible(False)
	
		plt.draw()
	
		ncols = np.round(len(lc_filenames)/4)
		legends = [lc_filename.split('_')[1].split('.')[0] for lc_filename in lc_filenames]
	        plt.legend(legends,loc=0,numpoints=1)
	        title_str = 'Light Curve ' + target_name
	#        plt.title(title_str)
		plt.rcParams['axes.labelsize'] = 20
		plt.rcParams['font.size'] = 20
	        plt.xlabel('JD (+'+str(t0)+')')
	        ylabel_str = 'magnitude'
	        plt.ylabel(ylabel_str)
	        plt.gca().invert_yaxis()
	#	ax = plt.gca()
	#	print ax.axis['left']
	#	ax.axis['left'].major_ticklabels.set_size(18)
	#	ax.axis['bottom'].major_ticklabels.set_size(18)
	#	ax2.axis['top'].major_ticklabels.set_size(14)
	
	#	ax.xaxis.get_label().set_fontproperties(font)
	#	ax.yaxis.get_label().set_fontproperties(font)
	#	plt.tick_params(labelsize= 18)
	#	ax.xaxis.set_tick_params(labelsize= 28)

		plt.grid()
	        plt.show()
	
	
	
