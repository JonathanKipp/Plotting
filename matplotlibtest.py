import matplotlib
import numpy as np
import post_proc
#print(matplotlib.__version__)
idx = 5
filename = "/Users/kipp/STB/output/weyl_theta_" + "{:d}/".format(idx) + "band_E.npy"
data = np.load(filename)
post_proc.plot_square_lattice(data)
