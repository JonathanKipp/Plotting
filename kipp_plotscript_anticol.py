#plots data from output file
import post_proc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm,colors
from matplotlib.backends.backend_pdf import PdfPages
from kipp_plotbib import canting_bandplotter
prefix = "/Users/kipp/STB/output/"
save_prefix = "/Users/kipp/STB/weyl_figures/Weyl_figures/"
filename = "path_rel_G-K-Kprime_outwards_anticol_anticol_theta"

canting_bandplotter(fname = filename,files = 20,bands = [0,1,2,3],combine = False,hall = True,save_str = 'anticol_inwards')
canting_bandplotter(fname = filename,files = 20,bands = [0,1,2,3],combine = True,hall = True,save_str = 'anticol_inwards')
