from kipp_plotbib import replacer, array_replacer,canting_bandplotter,replacer_samefile,sigma_trend_antisym
import numpy as np      
import glob, os, sys
from subprocess import call

def keyarrgen(rotation,checker,ext = 0.01,steps = 5,start = 0,start2 =0):
        print(rotation)
        if rotation == 'ferro':
                temp = np.pi/2*np.linspace(-1,2,4)
                temp2 = temp
                
        elif rotation == 'inplane':
                temp = np.pi/2*np.array([-1,1,-1,1])
                temp2 = np.pi/2*np.array([-1,1,1,-1])
                
        elif rotation == 'oop':
                temp = np.pi*np.array([0,1,0,1])
                temp2 = np.pi*np.array([0,1,1,0])
        elif rotation == 'scanoopferro':

                temp = np.pi*(0 + np.linspace(-ext,ext,steps))
                temp2 = - temp
        elif rotation == 'scanoopantiferro':

                temp = np.pi*(0 + np.linspace(-ext,ext,steps))
                temp2 = np.pi*(1 + np.linspace(-ext,ext,steps))
        elif rotation == 'scanipferro':

                temp = np.pi/2*(1 + np.linspace(0,ext,steps))
                temp2 = np.pi/2*(1 + np.linspace(0,ext,steps))
        elif rotation == 'scanipferroflip':

                temp = np.pi/2*(1 - np.linspace(0,ext,steps))
                temp2 = np.pi/2*(1 - np.linspace(0,ext,steps))                
        elif rotation == 'scanipantiferro1':

                temp = np.pi/2*(1 + np.linspace(-ext,ext,steps))
                temp2 = np.pi/2*(-1 - np.linspace(-ext,ext,steps))
        elif rotation == 'scanipantiferro2':
                temp = np.pi/2*(-1 - np.linspace(-ext,ext,steps))
                temp2 = np.pi/2*(1 + np.linspace(-ext,ext,steps))
        elif rotation == 'scanipantiferro3':
                temp = np.pi/2*(1 - np.linspace(-ext,ext,steps))
                temp2 = np.pi/2*(-1 + np.linspace(-ext,ext,steps))
        elif rotation == 'scanipantiferro4':
                temp = np.pi/2*(1 + np.linspace(-ext,ext,steps))
                temp2 = np.pi/2*(-1 - np.linspace(-ext,ext,steps))
        elif rotation == 'scancol':
                temp = np.pi*(start + np.linspace(0,ext,steps))
                temp2 = np.pi*(start2 + np.linspace(0,ext,steps))
        elif rotation == 'scancolflip':
                temp = np.pi*(-start + np.linspace(0,ext,steps))
                temp2 = np.pi*(-start2 + np.linspace(0,ext,steps))
        elif rotation == 'scansingular':
                temp = np.pi*(start + np.linspace(0,ext,steps))
                temp2 = np.pi*(start2 + np.zeros(steps))
        elif rotation == 'scansingularflip':
                temp = np.pi*(start - np.linspace(0,ext,steps))
                temp2 = np.pi*(start2 + np.zeros(steps))                
        else:
                print('Rotation not known. Valid rotations are: 0 ::"ferro", 1 :: "inplane", 2 :: "oop"')
                checker = 1
                return(None)
        returntemp = np.array([temp,temp2]).transpose()
        return(returntemp)
key = 'anticol_theta'
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
offset = 5./180
steps = 21
ext = 1.0
checker = 0
start1 = 1.0
start2 = 0.0
keyarr = keyarrgen(rotation,checker,ext = ext,steps = steps,start = start1 - offset,start2 = start2 + offset)
keyarr2 = keyarrgen(rotation2,checker,ext = ext,steps = steps,start = start1 - offset,start2 = start2 + offset)
fitarray = np.arange(steps/2-1,steps/2+2,dtype = int)
print(fitarray)

'''
        #REPLACE MAGTYPE
replacer_samefile(key = 'mag_type',keyarr = [magtype],prefix = prefix_inis,seedfname =seedfname)
        #REPLACE FERRO_PHI
replacer_samefile(key = 'ferro_phi',keyarr = [ferro_phi],prefix = prefix_inis,seedfname =seedfname)
        #REPLACE FERRO_THETA
replacer_samefile(key = 'ferro_theta',keyarr = [ferro_theta],prefix = prefix_inis,seedfname =seedfname)
        #REPLACE ANTICOL_PHI
replacer_samefile(key = 'anticol_phi',keyarr = [anticol_phi],prefix = prefix_inis,seedfname =seedfname)
        #REPLACE SAVE_STR
#replacer_samefile(key = 'band_prefix',keyarr = [save_prefix + ""],prefix = prefix_inis,seedfname =seedfname)
'''

if checker == 0:
        #REPLACE ANTICOL_THETA
        '''
        array_replacer(key = key,keyarr=keyarr,prefix = prefix_inis,seedfname = seedfname,savename = rotation)
        array_replacer(key = key,keyarr=keyarr2,prefix = prefix_inis,seedfname = seedfname,savename = rotation2)
        #RUN ALL SCRIPTS
        
        for file in glob.glob(prefix_inis + inputname + "*.cfg"):
	        try:
		        call (prefix_stb + "stb.x " + file,shell = 'True')
	
	        except:
		        print("Unexpected error:",sys.exc_info()[0])
        '''
        #PLOT DATA
        if plot_log == True:
                #canting_bandplotter(fname = inputname,fname2 = inputname2,files = keyarr.shape[0],bands = [0,1,2,3],combine = False,dos = False,hall = True,save_str = save_str,label = label)
                #canting_bandplotter(fname = inputname,files = keyarr.shape[0],bands = [0,1,2,3],combine = True,dos = True,hall = True,save_str = save_str,label = label)
                sigma_trend_antisym(fname = inputname,fname2 = inputname2,files = keyarr.shape[0],save_str = save_str,prom = prom,scale = 0.0,peakexp = 3,fitarray = fitarray,fitlog = False,coltype = 'col',hc_sym = 'antisym')
                sigma_trend_antisym(fname = inputname,fname2 = inputname2,files = keyarr.shape[0],save_str = save_str,prom = prom,scale = 0.0,peakexp = 3,fitarray = fitarray,fitlog = False,coltype = 'col',hc_sym = 'sym')