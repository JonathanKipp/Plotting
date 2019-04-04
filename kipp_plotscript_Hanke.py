import post_proc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm,colors
from matplotlib.backends.backend_pdf import PdfPages
from kipp_plotbib import bandplotter,E_reshaper

prefixHanke = '/Users/kipp/HankeCode'
save_prefixHanke = '/Users/kipp/HankeCode/output/'
code = 'tbcode'
fnameHanke = '/out_rashbabands'
save_str = 'path_rel_G-K-Kprime_Hanke'
fnameHankeberry = '/out_rashbaberryvvsum'
E,K = E_reshaper(fnamebands = fnameHanke,fnameberry = fnameHankeberry,nangels = 10,mode = 'line')

bandplotter(prefix = save_prefixHanke + code,fname = fnameHanke, bands = [0,1,2,3], save_str = save_str,files = 10,combine = True,hall = True)
bandplotter(prefix = save_prefixHanke + code,fname = fnameHanke, bands = [0,1,2,3], save_str = save_str,files = 10,combine = False,hall = True)

