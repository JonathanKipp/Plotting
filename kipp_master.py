from kipp_plotbib import replacer, array_replacer,canting_bandplotter,replacer_samefile,sigma_trend_antisym,keyarrgen
import numpy as np      
import glob, os, sys
from subprocess import call

def offset_in_deg(val):
        return(val/180.)
key = 'anticol_theta'
coltype = 'col'
rotation = 'scancol'
rotation2 = 'scancolflip'
seedfname = "path_rel_G-K-Kprime"
seedfnamesq = 'sq_fill_anticol'
save_str = seedfname + '_' + key + '_' + rotation
prefix_stb = "/Users/kipp/STB/"
prefix_inis = "/Users/kipp/STB/inis/"
save_prefix = "/output/" + seedfname
#theta = np.load(prefix + "m_theta.npy")
#keyarr = np.linspace(keymin,keymax,keystep)
band_prefix = ''
anticol_phi = [0.,0.]
ferro_theta = 0.0
ferro_phi = np.pi/2
magtype = 'anticol'
label = ["$\Gamma$", "$K^{'}$", "K", "$\Gamma$"]
plot_log = True
prom = 1.*10**-6

#inputname = "path_rel_G-K-Kprime_outwards_anticol_anticol"
inputname = seedfname + '_' + key + '_' + rotation
inputname2 = seedfname + '_' + key + '_' + rotation2

#inputname2 = 'sq_fill_anticol'

        #GENERATE REPLACEMENTS
offset = offset_in_deg(0)
steps = 4
ext = 0.5
checker = 0
start1 = 0.0
start2 = 0.0
keyarr = keyarrgen(rotation,checker,ext = ext,steps = steps,start1 = start1,start2 = start2 + offset)
keyarr2 = keyarrgen(rotation,checker,ext = ext,steps = steps,start1 = start1 + offset,start2 = start2)
lambdaarr = np.linspace(0,3,6)
tsoarr = np.linspace(0,3,6)
fitarray = np.arange(steps/2-1,steps/2+2,dtype = int)
if checker == 0:
        #REPLACE ANTICOL_THETA
        '''        #CHANGE ROTATION
        array_replacer(key = key,keyarr=keyarr,prefix = prefix_inis,seedfname = seedfname,savename = rotation)
        array_replacer(key = key,keyarr=keyarr2,prefix = prefix_inis,seedfname = seedfname,savename = rotation2)
                #CHANGE LAMBDA,TSO
        for j in range(keyarr.shape[0]):
                array_replacer(key = 'lambda',keyarr=lambdaarr,prefix = prefix_inis,seedfname = inputname + '_{:02}'.format(j),savename = '')
                array_replacer(key = 'lambda',keyarr=lambdaarr,prefix = prefix_inis,seedfname = inputname2 + '_{:02}'.format(j),savename = '')
                for l in range(lambdaarr.shape[0]):
                        array_replacer(key = 't_so',keyarr=tsoarr,prefix = prefix_inis,seedfname = inputname + '_{:02}'.format(j) + '_lambda' + '_{:02}'.format(l),savename = '')
                        array_replacer(key = 't_so',keyarr=tsoarr,prefix = prefix_inis,seedfname = inputname2 + '_{:02}'.format(j) + '_lambda' + '_{:02}'.format(l),savename = '')
        #RUN ALL SCRIPTS
        
        for file in glob.glob(prefix_inis + inputname + "*.cfg"):
	        try:
		        call ("mpirun -np 3 " + prefix_stb + "stb.x " + file,shell = 'True')
	
	        except:
		        print("Unexpected error:",sys.exc_info()[0])
        '''
        #PLOT DATA
        if plot_log == True:
                
                #canting_bandplotter(fname = inputname,fname2 = inputname2,files = keyarr.shape[0],bands = [0,1,2,3],combine = False,dos = True,hall = False,save_str = save_str,label = label,coltype = coltype)
                #canting_bandplotter(fname = inputname,files = keyarr.shape[0],bands = [0,1,2,3],combine = True,dos = True,hall = True,save_str = save_str,label = label)
                sigma_trend_antisym(fname = inputname,fname2 = inputname2,files = keyarr.shape[0],save_str = save_str,prom = prom,scale = 0.0,peakexp = 3,fitarray = fitarray,fitlog = False,coltype = coltype,hc_sym = 'antisym')
                sigma_trend_antisym(fname = inputname,fname2 = inputname2,files = keyarr.shape[0],save_str = save_str,prom = prom,scale = 0.0,peakexp = 3,fitarray = fitarray,fitlog = False,coltype = coltype,hc_sym = 'sym')            
        