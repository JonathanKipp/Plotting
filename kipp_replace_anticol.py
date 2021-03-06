#import numpy as np
#import re
from kipp_plotbib import replacer, array_replacer
import numpy as np
#atan = np.linspace(1.95, 2.05, 10)
#atan = np.arange(10,360,5)
#atan = np.logspace(-2, 1.0,12)
#rad = np.arange(21,60,3)

#n = 20				 # number of theta steps 
#theta = np.pi/2*np.flip(np.linspace(0,n),0)
# phi   = np.load("m_phi.npy")

#r = np.linspace(0.1, 0.9, 7)


#prefix = "/home/matthias/STB/inisv/"
prefix = "/Users/kipp/STB/inis/"
theta = np.load(prefix + "m_theta.npy")
seedfname = "path_rel_G-K-Kprime_outwards_anticol"
seedfnamesq = 'sq_fill_anticol'
cnt = 0
a = np.pi/2*np.linspace(0,2,20)
keyarr = np.array([-a,a]).transpose()

array_replacer(key ='anticol_theta',keyarr = keyarr,seedfname  = seedfname,prefix = prefix)
#array_replacer(key ='anticol_theta',keyarr = keyarr,seedfname  = seedfnamesq,prefix = prefix)

'''for num, val in enumerate(theta):
	with open(prefix + seedfname + '.cfg', "rt") as fin:
		with open(prefix + seedfname +"_{:02}.cfg".format(num), "wt") as fout:
			for line in fin:
                #line = line.replace('%m_phi%',   "{:011.6f}".format(phi[i]))
				line = line.replace('m_theta', "{:1.6f}".format(val))
				line = line.replace(seedfname,seedfname + "_{:02}".format(num))
                #line = line.replace('%middle%', "{:011.6f}".format(mid))
				fout.write(line)
'''