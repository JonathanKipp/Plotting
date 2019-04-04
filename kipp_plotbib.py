import post_proc
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm,colors,ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
import os
import re
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import scipy.constants as sp
from matplotlib.ticker import FormatStrFormatter
from numpy.polynomial.polynomial import polyval
from matplotlib.offsetbox import AnchoredText
prefix = "/Users/kipp/STB/output/"
save_prefix = "/Users/kipp/STB/weyl_figures/Weyl_figures/"
filename = "weyl_theta"
prefixHanke = '/Users/kipp/HankeCode'
save_prefixHanke = '/Users/kipp/HankeCode/output'
code = '/tbcode'
fnameHanke = '/out_rashbabands'
fnameHankeberry = 'out_rashbaberryvvsum'
nangles = 10
def format_func(value, tick_number):
    return("{:.2E}".format(value))
def replacer_samefile(key,keyarr,prefix,seedfname,save_str = 'w'):
	for num,val in enumerate(keyarr):
		with open(prefix + seedfname + '.cfg','r+') as fin:
			with open(prefix + seedfname + '_' + save_str + '.cfg','w+') as fout:
				for line in fin:
					pattern = key + "\s*=\s*['a-z]+|" + key + "\s*=\s*(\d+(.\d*))|" + key + "\s*=\s*"
					line = re.sub(pattern,re.sub('[\n\[\]]','',key + " = {}".format(val)),line)
					fout.write(line)
def linfunc(x,a,b,c):
    return(a*x**3 + b*x +c)
def replacer(key,keyarr,prefix,seedfname):
	#keyarr = np.linspace(keymin,keymax,keystep)
	for num,val in enumerate(keyarr):
		with open(prefix + seedfname + '.cfg','rt') as fin:
			with open(prefix + seedfname + '_' + key +"_{:02}.cfg".format(num), "wt") as fout:	
				for line in fin:
					#print(re.search('(?<=' + key + '\s=\s)',line))
					pattern = lambda:x + "\s*=\s*['a-z]+|" + x + "\s*=\s*(\d+(.\d*))|" + x + "\s*=\s*"
					line = re.sub(pattern,re.sub('[\n\[\]]','',key + " = {}".format(val)),line)
					line = re.sub(pattern)
					#line = re.sub('(?<=' + key + '\s=\s)',"{}".format(val),line)
					#line = re.sub('(?<=' + seedfname + ')','_' + key + "_{:02}".format(num),line)
					#line = line.replace(key, "{:1.6f}".format(val))
					#line = line.replace(seedfname,seedfname + "_{:02}".format(num))
					fout.write(line)
	return(seedfname + '_' + key)
def array_replacer(key,keyarr,prefix,seedfname,savename = ''):
	#keyarr = np.linspace(keymin,keymax,keystep)
	#keyarr = np.array(x + step for step in )
	keyarrstr = [re.sub('[\n\[\]]','',np.array2string(x,precision = 4, separator = ',')) for x in keyarr]
	for num,val in enumerate(keyarrstr):
		with open(prefix + seedfname + '.cfg','rt') as fin:
			with open(prefix + seedfname + '_' + key + '_' + savename +"_{:02}.cfg".format(num), "wt") as fout:	
				for line in fin:
					#print(re.search('(?<=' + key + '\s=\s)',line))
					line = re.sub('(?<=' + key + '\s=\s)',val,line)
					#line = re.sub('(?<=' + key + '\s=\s)',"{:1.6f}".format(val),line)
					line = re.sub('(?<=' + seedfname + ')','_' + key + '_' +savename + "_{:02}".format(num),line)
					#line = line.replace(key, "{:1.6f}".format(val))
					#line = line.replace(seedfname,seedfname + "_{:02}".format(num))
					fout.write(line)
				fout.close()
			
	


def versionmaker(prefix,fname):
	today = date.isoformat(date.today())
	version = 0
	if os.path.isfile(prefix + fname + 'today.npy')==True:
		today_old  = np.load(prefix + fname + 'today.npy')
		if today != today_old[0]:
			np.save(prefix + fname + 'today',[today])
		elif today==today_old:
			if os.path.isfile(prefix + fname + 'version.npy')==True:
				version = np.load(prefix + fname + 'version.npy')[0] + 1
				np.save(prefix + fname + 'version.npy',[version])
			else:
				np.save(prefix + fname + 'version.npy',[version])
	else: 
		np.save(prefix + fname + 'today.npy',[today])
		np.save(prefix + fname + 'version.npy',[0])
	return(today, version)
def E_reshaper(prefix = prefixHanke,save_prefix = save_prefixHanke, code = code,nangels = nangles,fnamebands = fnameHanke,fnameberry = fnameHankeberry,sfname = '',mode = 'square'):
	rbands = np.loadtxt(prefix + code + fnamebands)
	berry = np.loadtxt(prefix + code + fnameberry)
	if mode == 'square':
		sq_shape = int(np.sqrt(rbands.shape[0]/4/nangles))
		E = rbands[:,4].reshape(nangles,4,sq_shape,sq_shape)
		K = rbands[:,1:3].reshape(10,4,sq_shape,sq_shape,2)
		#E_sig = berry[:,0].reshape(nangles,sq_shape,sq_shape)
		#sigma = berry[:,2:].reshape(nangles,sq_shape,sq_shape,4)
	elif mode == 'line':
		sq_shape = int(rbands.shape[0]/4/nangles)
		berry_shape = int(berry.shape[0]/nangles)
		E = rbands[:,3].reshape(nangles,4,sq_shape)
		K = rbands[:,1:3].reshape(nangles,4,sq_shape,2)
		E_sig = berry[:,0].reshape(nangles,berry_shape)
		sigma = berry[:,2:].reshape(nangles,berry_shape,4)
	for j in range(nangels):
		np.save(save_prefix + code + fnamebands + sfname + '_{:02}/'.format(j) + 'band_E',E[j])
		np.save(save_prefix + code + fnamebands + sfname + '_{:02}/'.format(j) + 'band_k',K[0,0])
		np.save(save_prefix + code + fnamebands + sfname + '_{:02}/'.format(j) + 'hall_cond_E',E_sig[j])
		np.save(save_prefix + code + fnamebands + sfname + '_{:02}/'.format(j) + 'hall_cond',sigma[j,:,1])
	return(E,K)

def bandplotter(fname,files,prefix = '/Users/kipp/STB/output/',saveprefix = '/Users/kipp/STB/weyl_figures/',save_str = None,combine =True,bands = [1,2],max_idx = 10, hall = False): #combine = False: plot each bandstructure in one figure
	norm = colors.Normalize(vmin = 0, vmax = files)
	n_kpts = 300
	ticks =post_proc.calc_tick_pos(3,n_kpts)
	label = ["$\Gamma$", "$K^{'}$", "K", "$\Gamma$"]
	if combine == True:
		d,v = versionmaker(prefix,fname + 'combine')
		fig,(axis,axis2,axes3) = plt.subplots(nrows = 1,ncols = 3,constrained_layout = True)
		axis.set_xticks(ticks)
		axis.set_xticklabels(label)
		axis.set_ylabel('Energy [eV]')
		axis.set_title('Bandstructure for $\Theta$ in $[0,\pi/2]$')
		axis.axvline(ticks[1],color = 'k',linewidth = '0.5')
		axis.axvline(ticks[2],color = 'k',linewidth = '0.5')
		axis2.set_ylabel('Energy [eV]')
		axis2.set_xlabel('$\sigma_{xy}$')
		axis2.set_title('Hall conductance')
		for j in range(files):
			E = np.load(prefix + fname + '_{:02}/'.format(j) + 'band_E.npy')[bands]
			k = np.arange(0,E.shape[1])
			if hall == True:
				hc_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond_E.npy')
				hc = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond.npy')
				axis2.plot(hc,hc_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(j)))
			#else:
			#	axis2.remove()
			for e in E:
				axis.plot(k,e,lw=0.5,color = cm.get_cmap('plasma')(norm(j)))
		sm = plt.cm.ScalarMappable(cmap = 'plasma',norm = norm)
		sm._A = []
		cbar = fig.colorbar(cax = axis2,ax = [axes,axes2,axes3],mappable = sm,ticks = [0,files])
		cbar.set_ticklabels(["0","$\pi$"])
		cbar.set_label("$\Theta_{ferro}$")

		if save_str == None:
			with PdfPages(save_prefix + d + '_' + str(v) + '_' + fname + 'bands.pdf') as pdf:
				pdf.savefig(fig)
		else:
			with PdfPages(save_prefix + d + '_' + str(v) + '_' + save_str + 'bands.pdf') as pdf:
				pdf.savefig(fig)
	if combine == False:
		d,v = versionmaker(prefix,fname + 'notcombine')
		for j in range(files):
			fig, (axis,axis2) = plt.subplots(nrows = 1, ncols = 2,constrained_layout = True)
			axis.set_xticks(ticks)
			axis.set_xticklabels(label)
			axis.set_ylabel('Energy [eV]')
			axis.set_title('Bandstructure for $\Theta = {:1.2f}$ degrees'.format(j/(files-1)*180))
			axis.axvline(ticks[1],color = 'k',linewidth = '0.5')
			axis.axvline(ticks[2],color = 'k',linewidth = '0.5')
			axis2.set_ylabel('Energy [eV]')
			axis2.set_xlabel('$\sigma_{xy}$')
			axis2.set_title('Hall conductance')
			E = np.load(prefix + fname + '_{:02}/'.format(j) + 'band_E.npy')[bands]
			k = np.arange(0,E.shape[1])
			for e in E:
				axis.plot(k,e,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
			if hall == True:
				hc_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond_E.npy')
				hc = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond.npy')
				axis2.plot(hc,hc_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
			else:
				axis2.remove()
			if save_str == None:
				with PdfPages(save_prefix + d + '_' + str(v) + '_' + fname + '_{:02}.pdf'.format(j)) as pdf:
					pdf.savefig(fig)
			else:
				with PdfPages(save_prefix + d + '_' + str(v) + '_' + save_str + '_{:02}.pdf'.format(j)) as pdf:
					pdf.savefig(fig)
			plt.close(fig)

def canting_bandplotter(fname,fname2,files,prefix = '/Users/kipp/STB/output/',saveprefix = '/Users/kipp/STB/weyl_figures/',save_str = None,combine =True,bands = [1,2],max_idx = 10, dos = False,hall = False,cbarlabels = [0,np.pi],label = ["$\Gamma$", "$K^{'}$", "K", "$\Gamma$"]): #combine = False: plot each bandstructure in one figure
	norm = colors.Normalize(vmin = 0, vmax = files)
	n_kpts = 300
	ticks =post_proc.calc_tick_pos(3,n_kpts)
	#label = ["$\Gamma$", "$K^{'}$", "K", "$\Gamma$"]
	if combine == True:
		d,v = versionmaker(prefix,fname + 'combine')
		m_phi = np.load(prefix + fname + '_{:02}/'.format(0) + 'm_phi.npy')
		thetas = []
		numplotsx = 3
		numplotsy = 2
		xsize = 6.4*numplotsx
		ysize = 4.8*numplotsy
		sizetup = (xsize,ysize)
		fig,ax = plt.subplots(nrows = numplotsy,ncols = numplotsx,sharey = 'row',figsize = sizetup,constrained_layout = True)
		#axis.set_aspect(100)
		#axis2.set_aspect(1./4)
		#axis3.set_aspect(1./8)
		ax[0,0].set_xticks(ticks)
		ax[0,0].set_xticklabels(label)
		ax[0,0].set_ylabel('Energy [eV]')
		ax[0,0].set_title('Bandstructure')
		ax[0,0].axvline(ticks[1],color = 'k',linewidth = '0.5')
		ax[0,0].axvline(ticks[2],color = 'k',linewidth = '0.5')
		ax[0,1].set_ylabel('Energy [eV]')
		ax[0,1].set_xlabel('$\sigma_{xy}$')
		ax[0,1].set_title('Hall conductance')
		ax[0,2].set_title('DOS')
		ax[1,0].set_xticks([])
		ax[1,0].set_yticks([])
		textstr = "Input Parameters:\n" + "$\Phi_{cant,1}$" + " = {0:1.2f}".format(m_phi[0]) + "\n$\Phi_{cant,2}$" + " = {0:1.2f}\n".format(m_phi[1]) + "Path = {0} - {1} - {2} - {3}".format(*label)
		ax[1,0].text(0, 0, textstr,size = 18,
			horizontalalignment='left',
			verticalalignment='bottom')
			
		ax[1,0].set_axis_off()
		ax[1,-2].remove()
		ax[1,-1].remove()		
		axins = inset_axes(ax[1,0], width=1.3, height=0.9,loc = 'upper left')
		axins.tick_params(tick1On = False,labelleft=False, labelbottom=False)
		axins.set_title('Spin Configuration')
		axins.set_xlim(-1.5,2.5)
		axins.set_ylim(-1.2,1.2)
		for j in range(files):
			m_theta = np.load(prefix + fname + '_{:02}/'.format(j) + 'm_theta.npy')
			thetas = np.append(thetas,m_theta[0])
			dx = 0.5*np.sin(m_theta)
			dz = 0.5*np.cos(m_theta)
			axins.arrow(1.15, 0, dx[0],dz[0], head_width=0.1, head_length=0.1,color = cm.get_cmap('plasma')(norm(j)))
			axins.arrow(-0.15, 0, dx[1],dz[1], head_width=0.1, head_length=0.1,color = cm.get_cmap('plasma')(norm(j)))	
			E = np.load(prefix + fname + '_{:02}/'.format(j) + 'band_E.npy')[bands]
			DOS = np.load(prefix + fname + '_{:02}/'.format(j) + 'DOS.npy')
			DOS_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'DOS_E.npy')
			k = np.arange(0,E.shape[1])
			#else:
			#	axis2.remove()
			for e in E:
				ax[0,0].plot(k,e,lw=0.5,color = cm.get_cmap('plasma')(norm(j)))
			if hall == True:
				hc_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond_E.npy')
				hc = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond.npy')
				ax[0,1].plot(hc,hc_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(j)))
			else:
				ax[0,1].remove()
			ax[0,2].plot(DOS,DOS_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(j)))
		sm = plt.cm.ScalarMappable(cmap = 'plasma',norm = norm)
		sm._A = []
		cbar = fig.colorbar(ax = ax[1,0],mappable = sm,ticks = [0,files],aspect = 40)
		cbar.set_ticklabels(['{:1.2f}'.format(x) for x in [min(thetas),max(thetas)]])
		cbar.set_label("$\Theta_{anticol}$")
		if save_str == None:
			with PdfPages(save_prefix + d + '_' + str(v) + '_' + fname + 'bands.pdf') as pdf:
				pdf.savefig(fig)
		else:
			with PdfPages(save_prefix + d + '_' + str(v) + '_' + save_str + 'bands.pdf') as pdf:
				pdf.savefig(fig)
	if combine == False:
		d,v = versionmaker(prefix,fname + 'notcombine')
		numplotsx = 4
		numplotsy = 2
		xsize = 6.4*numplotsx
		ysize = 4.8*numplotsy
		sizetup = (xsize,ysize)
		for j in range(files):
			m_theta = np.load(prefix + fname + '_{:02}/'.format(j) + 'm_theta.npy')
			m_phi = np.load(prefix + fname + '_{:02}/'.format(j) + 'm_phi.npy')
			key = 'temperature'
			pattern = key + "\s*=\s*['a-z]+|" + key + "\s*=\s*(\d+(.\d*))|" + key + "\s*=\s*"
			with open(prefix + fname + '_{:02}/'.format(j) + 'setup.cfg','r+') as fin:
				for line in fin:
					s = re.search(pattern,line)
					if s!=None:
						if s.group(2)!=None:
							temperature = float(s.group(2))*1000
			#temperature = np.load(prefix + fname + '_{:02}/'.format(j) + 'temperature.npy')                    
			fig, ax = plt.subplots(nrows = numplotsy, ncols = numplotsx,sharey = 'row',figsize = sizetup,constrained_layout = True)
			ax[0,0].set_xticks(ticks)
			ax[0,0].set_xticklabels(label)
			ax[0,0].set_ylabel('Energy [eV]')
			ax[0,0].set_title('Bandstructure')
			ax[0,1].set_title('Hall conductance')
			ax[0,2].set_title('Hall conductance antisymmetric')
			ax[0,3].set_title('Hall conductance symmetric')
			ax[0,0].axvline(ticks[1],color = 'k',linewidth = '0.5')
			ax[0,0].axvline(ticks[2],color = 'k',linewidth = '0.5')
			ax[0,1].set_xlabel('$\sigma_{xy}$')
			ax[0,2].set_xlabel('$\sigma_{xy,antisym}$')
			ax[0,3].set_xlabel('$\sigma_{xy,sym}$')
			ax[1,0].set_xticks([])
			ax[1,0].set_yticks([])
			ax[1,0].set_axis_off()
			thes = np.pi/2*np.arange(-1,3,1)
			thesdiff1 = thes - m_theta[0]
			thesdiff2 = thes - m_theta[1]
			thesmin1 = 180/np.pi*thes[np.argmin(abs(thesdiff1))]
			thesmin2 = 180/np.pi*thes[np.argmin(abs(thesdiff2))]
			minthesdiff1 = 180/np.pi*thesdiff1[np.argmin(abs(thesdiff1))]
			minthesdiff2 = 180/np.pi*thesdiff2[np.argmin(abs(thesdiff2))]
			#"Temperature = {:1.2f}".format(temperature)            
			textstr = "Input Parameters:\nTemperature = {:1.2f} K\n".format(temperature) + "$\Theta_{cant,1}$" + " = {0:1.2f} + {1:1.2f} deg.".format(thesmin1,minthesdiff1) + "\n$\Theta_{cant,2}$" + " = {0:1.2f} + {1:1.2f} deg.".format(thesmin2,minthesdiff2) + "\nPath = {0} - {1} - {2} - {3}".format(*label)
			ax[1,0].text(0, 0, textstr,size = 18,
				horizontalalignment='left',
				verticalalignment='bottom')
			ax[1,-1].remove()
			ax[1,-2].remove()
			
			axins = inset_axes(ax[1,0], width=1.3, height=0.9,loc = 'upper left')
			axins.tick_params(tick1On = False,labelleft=False, labelbottom=False)
			dx = 0.5*np.sin(m_theta)
			dz = 0.5*np.cos(m_theta)
			axins.arrow(1.2, 0, dx[0],dz[0], head_width=0.1, head_length=0.1)
			axins.arrow(-0.2, 0, dx[1],dz[1], head_width=0.1, head_length=0.1)
			axins.set_xlim(-1.2,2.2)
			axins.set_ylim(-1.2,1.2)
			axins.set_title('Spin Configuration')
			E = np.load(prefix + fname + '_{:02}/'.format(j) + 'band_E.npy')[bands]
			k = np.arange(0,E.shape[1])
			if dos == True:
				DOS = np.load(prefix + fname + '_{:02}/'.format(j) + 'DOS.npy')
				DOS_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'DOS_E.npy')
			for e in E:
				ax[0,0].plot(k,e,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
			if hall == True:
				hc_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond_E.npy')
				hc = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond.npy')
				hc_2 = np.load(prefix + fname2 + '_{:02}/'.format(j) + 'hall_cond.npy')
				for l in range(1,4):
					ax[0,l].xaxis.set_major_locator(ticker.AutoLocator())
					ax[0,l].xaxis.set_major_formatter(ticker.ScalarFormatter())
					ax[0,l].xaxis.get_major_formatter().set_scientific(True)
					ax[0,l].xaxis.get_major_formatter().set_powerlimits(lims = (-3,4))
				ax[0,1].plot(hc,hc_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
				ax[0,2].plot((hc-hc_2)/2,hc_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
				ax[0,3].plot((hc+hc_2)/2,hc_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
				#axis2.set_yticks(axis.get_yticks())
			else:
				ax[0,1].remove()
				ax[0,2].remove()
				ax[1,3].remove()
			if dos == True:
				ax[1,1].plot(DOS,DOS_E,lw = 0.5,color = cm.get_cmap('plasma')(norm(0)))
			else:
				ax[1,1].remove()
			if save_str == None:
				with PdfPages(save_prefix + d + '_' + str(v) + '_' + fname + '_{:02}.pdf'.format(j)) as pdf:
					pdf.savefig(fig)
			else:
				with PdfPages(save_prefix + d + '_' + str(v) + '_' + save_str + '_{:02}.pdf'.format(j)) as pdf:
					pdf.savefig(fig)
			plt.close(fig)

def sigma_trend_antisym(fname,fname2,files,prefix = '/Users/kipp/STB/output/',saveprefix = '/Users/kipp/STB/weyl_figures/',save_str = None,prom = 0.03,scale = 0.0, fitlog = False,coltype = 'noncol',hc_sym = 'antisym',peakexp = 0,fitarray = []):
	myfsize = 26
	mymarkersize = 12
	E_points = 3000
	width = 100
	widthuse = width
	promuse = prom
	#E_default = np.arange(1,E_points,E_points/4,dtype = int)
	E_default = np.array([580, 1060, 1530],dtype = int)
	d,v = versionmaker(prefix,fname + 'sigmatrend')
	norm = colors.Normalize(vmin = 0, vmax = files)
	peakslist = np.zeros((files,peakexp + 2))
	peaksposlist = np.tile(E_default,(files,1))
	peak_E = np.zeros((files,peakexp))
	hall_quant = (sp.e)**2/sp.h
	fitlogarr = np.zeros(files,dtype = int)
	fitlogarr[fitarray] = 1
	key = 'temperature'
	pattern = key + "\s*=\s*['a-z]+|" + key + "\s*=\s*(\d+(.\d*))|" + key + "\s*=\s*"
	with open(prefix + fname + '_{:02}/'.format(0) + 'setup.cfg','r+') as fin:
		for line in fin:
			s = re.search(pattern,line)
			if s!=None:
				if s.group(2)!=None:
					temperature = float(s.group(2))*1000
	for j in range(files):
		hc_E = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond_E.npy')
		hc = np.load(prefix + fname + '_{:02}/'.format(j) + 'hall_cond.npy')
		hc_opp = np.load(prefix + fname2 + '_{:02}/'.format(j) + 'hall_cond.npy')
		m_theta = 180./np.pi*np.load(prefix + fname + '_{:02}/'.format(j) + 'm_theta.npy')
		if hc_sym == 'sym':
			hc_fin = (hc + hc_opp)/2
		elif hc_sym == 'antisym':
			hc_fin = (hc - hc_opp)/2
			widthuse = 1
			promuse = 0.01
		elif hc_sym == 'bare':
			hc_fin = hc 
		#temperature = np.load(prefix + fname + '_{:02}/'.format(j) + 'temperature.npy')        
		#bands = np.load(prefix + inputname + '_{:02}/'.format(j) + 'band_E.npy')
		if coltype == 'noncol':        
			theta_int = m_theta[0]- m_theta[1]
		elif coltype == 'col':
				theta_int = (m_theta[0] + m_theta[1])/2
		length = int(hc_fin.shape[0]/2+width+10)
		peaksinput = abs(hc_fin[:length])
		peaks, props = find_peaks(peaksinput, prominence = (1+scale)*prom,width = width,distance = width)
		#print(peaks,props['prominences'],props['widths'])
		peakslist[j,1] = theta_int
		if peaks.shape == peakslist[j,2::].shape:       
			#peakslist[j,0] = theta_int
			peaksposlist[j,:] = peaks
			peakslist[j,2::] = hc_fin[peaks]
			peakslist[j,0] = 1
			peak_E[j] = hc_E[peaks]
		else:
			peak_E[j] = hc_E[E_default]
			peakslist[j,2::] = hc_fin[E_default]
			peakslist[j,0] = 1
			print("j = {}\nprom = {}\ntheta_int = {}\n".format(j,(1+scale)*prom,theta_int))
			print(peaks)#,hc_fin[peaks],hc_E[peaks])
	#print(peakslist)
	angles = peakslist[:,1]
	if all(peakslist[:,0] == 0.0) == False:
		fitin = peakslist[fitlogarr==1]        
		peakslist = peakslist[peakslist[:,0]!=0.0]
		fitin = fitin[fitin[:,0]!=0.0]
		#print(peakslist.shape)
	else:
		fitin = peakslist[fitlogarr==1]       
	print(peakslist.shape)
	npeaks = peakslist.shape[1]-2
	print(npeaks)
	nyplots = 2
	xsize = 2*6.4*npeaks
	ysize = 2*4.8*nyplots
	sizetup = (xsize,ysize)
	fig,ax = plt.subplots(nrows = nyplots,ncols = npeaks,figsize = sizetup,constrained_layout = True)
	#ax[1,-1].remove()
	for l in range(npeaks):       
		#fig,ax = plt.subplots(nrows = 1,ncols = 1,figsize = sizetup,constrained_layout = True)        
		if coltype == 'noncol':       
			labelstr = '$\sigma_{asym,fit} = a (\Theta_{1}-\Theta_{2}) + b$\n' 
			labelstr2 = '$\sigma_{asym,noncol}(\Theta{1}-\Theta{2},E)$'
			labelstr3 = r"$\frac{\Delta\sigma_{asym,noncol}}{\Delta\Theta}(\Theta{1}-\Theta{2},E)$"
		elif coltype == 'col':
			labelstr2 = '$\sigma_{asym,col}(\Theta_{col},E)$'
			labelstr3 = r"$\frac{\Delta\sigma_{asym,col}}{\Delta\Theta}(\Theta_{col},E)$"                
			labelstr = '$\sigma_{asym,fit} = a (\Theta_{col}) + b$\n'
		if fitlog == True:
			#fitinx = peakslist[fitlogarr==1,1]
			#fitiny = peakslist[fitlogarr==1,l+2]            
			p = np.polyfit(fitin[:,1],fitin[:,l+2],deg = 1) 
			labelstr = labelstr + r'a = {:1.2f}'.format(p[0])+ r'[$\frac{e^2}{h}$ per deg.]'+'\nb = {:1.2f}'.format(p[1]) + r'[$\frac{e^2}{h}$]'
		delta = 10**(-10)
		hc_pwd = (peakslist[:,l+2]-np.roll(peakslist[:,l+2],1))/(peakslist[:,1]-np.roll(peakslist[:,1],1) + delta)
		hc_pwd2 = (peakslist[:,l+2]-np.roll(peakslist[:,l+2],-1))/(peakslist[:,1]-np.roll(peakslist[:,1],-1) + delta)
		hc_pwd = (hc_pwd+hc_pwd2)/2
		hc_pwd[[0,-1]] = hc_pwd[[1,-2]]
		ax[0,l].grid(True)
		dist = max(peakslist[:,l+2])-min(peakslist[:,l+2])
		step = max(peakslist[:,l+2])/2
		nsteps = 4
		exp = 0
		if dist<1:
			for power in range(10):
				y = dist*10**power
				if y>=1:
					exp = power
					ndist = np.round(dist*10**(exp)/nsteps,1)
					step = ndist*10**(-exp)
					break
		ax[0,l].yaxis.set_major_locator(plt.MultipleLocator(step))
		ax[0,l].yaxis.set_major_formatter(plt.FuncFormatter(format_func))
		xstep = 30
		ax[0,l].xaxis.set_major_locator(plt.MultipleLocator(xstep))
		ax[0,l].xaxis.set_major_formatter(plt.ScalarFormatter())
		#ax[0,l].yaxis.set_major_locator(ticker.AutoLocator())
		#ax[0,l].yaxis.set_major_formatter(ticker.ScalarFormatter())
		#ax[0,l].yaxis.get_major_formatter().set_scientific(True)
		#ax[0,l].yaxis.get_major_formatter().set_powerlimits(lims = (-3,4))        
		ax[0,l].tick_params(labelsize = myfsize)
		ax[0,l].plot(peakslist[:,1],peakslist[:,l+2],'b',markersize = mymarkersize,label = labelstr2)
		ax[0,l].plot(peakslist[:,1],hc_pwd,'g',markersize = mymarkersize,label = labelstr3)
		if fitlog == True:
			yfit = polyval(fitin[:,1],np.flip(p))
			ax[0,l].plot(fitin[:,1],yfit,'r',linewidth = 1.,label = labelstr)
		ax[0,l].set_title('Peak at $E = {:1.02f}$ eV'.format(peak_E[0,l]),fontsize = myfsize)
		if coltype == 'noncol':       
			ax[0,l].set_xlabel('$\Theta_{1}-\Theta_{2}$ [deg.]',fontsize = myfsize)
		elif coltype == 'col':        
			ax[0,l].set_xlabel('$\Theta_{col}$ [deg.]',fontsize = myfsize)
		ax[0,l].set_ylabel(r'$\sigma_{asym,xy}$ [$\frac{e^2}{h}$]',fontsize = myfsize)
		ax[0,l].legend(loc = 'best',fontsize = myfsize,frameon = False)
		#if len(peakslist[:,0])<10:
		#	ax[0,l].set_xticks(peakslist[:,1])
		#ax[0,l].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	for j in [0,1,2]:
		ax[1,j].tick_params(labelsize = myfsize-2)
		ax[1,j].grid(True)
		index = int((files-1)/2*j)
		hc_E = np.load(prefix + fname + '_{:02}/'.format(index) + 'hall_cond_E.npy')
		hc = np.load(prefix + fname + '_{:02}/'.format(index) + 'hall_cond.npy')
		hc_opp = np.load(prefix + fname2 + '_{:02}/'.format(index) + 'hall_cond.npy')
		if hc_sym == 'sym':
			symlabel = r'$\sigma_{sym,xy}$ [$\frac{e^2}{h}$]'    
			symtitle = 'Symmetric part $\sigma_{sym} = (\sigma(\Theta) + \sigma(-\Theta))/2$'
			symstr = '$(\Theta_{1} + \Theta_{2})/2 = $'
			hc_fin = (hc + hc_opp)/2
		elif hc_sym == 'antisym':
			hc_fin = (hc - hc_opp)/2
			symlabel = r'$\sigma_{asym,xy}$ [$\frac{e^2}{h}$]'    
			symtitle = 'Antisymmetric part $\sigma_{asym} = (\sigma(\Theta) - \sigma(-\Theta))/2$'            
			symstr = '$(\Theta_{1} - \Theta_{2})/2 = $'            
		elif hc_sym == 'bare':
			hc_fin = hc
			symlabel = r'$\sigma_{xy}$ [$\frac{e^2}{h}$]'    
			symtitle = 'Bare part $\sigma = (\sigma(\Theta)$'            
			symstr = '$\Theta = $'            
		textstr = 'Temperature = {:1.2f} K\n'.format(temperature) + symstr + '{:1.2f} deg.'.format(angles[index])
		dist = max(hc_fin)-min(hc_fin)
		step = 1.0
		nsteps = 4
		exp = 0
		if dist<1:
			for power in range(10):
				y = dist*10**power
				if y>=1:
					exp = power
					ndist = np.round(dist*10**(exp)/nsteps,1)
					step = ndist*10**(-exp)
					break
		ax[1,j].xaxis.set_major_locator(plt.MultipleLocator(step))
		ax[1,j].xaxis.set_major_formatter(plt.FuncFormatter(format_func))
		#ax[1,j].xaxis.get_major_formatter().set_scientific(True)
		#ax[1,j].xaxis.get_major_formatter().set_powerlimits(lims = (-3,4))
		ax[1,j].plot(hc_fin,hc_E,color = 'blue',label = textstr)
		ax[1,j].legend(loc = 'lower left',fontsize = myfsize,frameon = False, edgecolor = 'k')
		ax[1,j].set_ylabel('E [eV]',fontsize = myfsize)
		ax[1,j].set_xlabel(symlabel,fontsize = myfsize)
		ax[1,j].set_title(symtitle,fontsize = myfsize)
		#ax[1,j].set_xticks(max(abs(hc_fin))*np.arange(-1,1.1,0.5))
		for p in peaksposlist[index]:
			#textpos = (-np.sign(hc[p])*0.2,hc_E[p])
			textpos = (0.0,hc_E[p])            
			ax[1,j].axhline(hc_E[p],color = 'k',linewidth = '0.5')
			an = ax[1,j].annotate("{:1.2f} eV".format(hc_E[p]), xy=(0., 0.),
			xytext=textpos, textcoords="data",
			va="bottom", ha="left",
			bbox=dict(boxstyle="round", fc="w"),fontsize = myfsize
			)
		#if hc_sym == 'antisym':            
		#	axins = inset_axes(ax[1,j], width=2.6, height=1.8,loc = 'center right')
		#else:
		#axins = inset_axes(ax[1,j], width=2.6, height=1.8,loc = 'lower right')            
		axins = inset_axes(ax[1,j], width=2.6, height=1.8,loc = 'lower right')
		axins.tick_params(tick1On = False,labelleft=False, labelbottom=False)
		axins.set_title('Spin Configuration',fontsize = myfsize)
		axins.set_xlim(-1.8,2.8)
		axins.set_ylim(-1.2,1.2)
		m_theta = np.load(prefix + fname + '_{:02}/'.format(int((files-1)/2*j)) + 'm_theta.npy')
		m_theta_opp = np.load(prefix + fname2 + '_{:02}/'.format(int((files-1)/2*j)) + 'm_theta.npy')
		dx = np.sin(m_theta)
		dz = np.cos(m_theta)
		dx_opp = np.sin(m_theta_opp)
		dz_opp = np.cos(m_theta_opp)        
		#print(m_theta,dx,dz)
		axins.arrow(1.65, 0, dx[0],dz[0], head_width=0.1, head_length=0.1,color = 'blue')
		axins.arrow(-0.65, 0, dx[1],dz[1], head_width=0.1, head_length=0.1,color = 'blue')
		axins.arrow(1.65, 0, dx_opp[0],dz_opp[0], head_width=0.1, head_length=0.1,color = 'red')
		axins.arrow(-0.65, 0, dx_opp[1],dz_opp[1], head_width=0.1, head_length=0.1,color = 'red')
	if save_str == None:
		with PdfPages(save_prefix + d + '_' + str(v) + '_' + fname + '_' + coltype + '_' + hc_sym + '_sigmatrend.pdf') as pdf:
			pdf.savefig(fig)
	else:
		with PdfPages(save_prefix + d + '_' + str(v) + '_' + save_str + '_' + coltype + '_' + hc_sym + '_sigmatrend.pdf') as pdf:
			pdf.savefig(fig)
	plt.close(fig)

