import numpy as np

prefix = 'Users/kipp/HankeCode'
save_prefix = 'Users/kipp/HankeCode/output'
code = '/tbcode'
fname = '/out_rashbabands'
def E_reshaper(prefix = prefix,save_prefix = save_prefix, code = code,fname,sfname):
    rbands = np.loadtxt(prefix + code + fname)
    sq_shape = int(np.sqrt(rbands.shape[0]/4))
    E = rbands[:,4].reshape(4,sq_shape,sq_shape)
    K = rbands[:,1:3].reshape(sq_shape,sq_shape,2)
    np.save(save_prefix + sfname + '/band_E',E)
    np.save(save_prefix + sfname + '/band_k',K)