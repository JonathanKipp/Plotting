import matplotlib.gridspec as gridspec
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import find_peaks,argrelextrema
from datetime import date
from numpy.polynomial.polynomial import polyval
from matplotlib import cm,colors,ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

def reshape_func(sq1,sq2,E,E_shape):
    
    #       RESHAPE 1D-ARRAY INTO SQ1xSQ2 2D-ARRAY AND RETURN ENTRIES ALONG THE 
    #       PATH ALONG THE INDICES SPECIFIED

    x = np.arange(230,1570                             ,dtype = 'int')
    y = 1522*np.ones(x.shape[0]                        ,dtype = 'int')
    k1 = np.arange(x.shape[0]                          ,dtype = 'int')
    x2 = np.arange(1569,901,-1                         ,dtype = 'int') 
    y2 = np.linspace(1522,3045,x2.shape[0]             ,dtype = 'int')
    k2 = k1[-1] + 1 + np.arange(x2.shape[0]            ,dtype = 'int')
    x3 = np.arange(901,236,-1                          ,dtype = 'int') 
    y3 = np.linspace(3045,1522,x3.shape[0]             ,dtype = 'int')
    k3 = k2[-1] + 1 + np.arange(x3.shape[0]            ,dtype = 'int')
    x = np.append(x,[x2])
    y = np.append(y,[y2])
    x = np.append(x,[x3])
    y = np.append(y,[y3])
    E_sq = E.reshape(E_shape,sq1,sq2)
    E_line = E_sq[:,y,x]
    return(E_line)

def format_func(value, tick_number):
    
    #       SWITCH TICK FORMAT TO POWER FORMAT IF SMALLER THAN THRESH
    
    thresh = 10**(-4)
    if abs(value) <= thresh and abs(value)!= 0.0:
        return("{:.2E}".format(value))
    else:
        return("{:1.4f}".format(value))

def hc_fin_func(hc,hc_opp,hc_sym):
    
    #       RETURN SYMMETRIC OR ANTISYMMETRIC PART OF HC

    if hc_sym == 'sym':
        hc_fin = (hc + hc_opp)/2
    elif hc_sym == 'antisym':
        hc_fin = (hc - hc_opp)/2
    elif hc_sym == 'bare':
        hc_fin = hc
    return(hc_fin)

def col_type_func(m_theta,case):
    
    #       CASE == 0:                     THETA1=-THETA2 = THETA_NC
    #       CASE == 1:              THETA1= 0, THETA2 = +/- THETA_NC
    #       CASE == 2:       THETA1= 0/THETA_NC, THETA2 = THETA_NC/0
    #       CASE == 3:     THETA1 = +/-THETA_NC THETA2 = -/+THETA_NC
    #       IN DEGREES
    
    if case == 0:
        theta_nc =             (m_theta[0] - m_theta[1])/2
        theta_col =            (m_theta[0] + m_theta[1])/2
    elif case == 1:
        theta_nc =                 m_theta[0] - m_theta[1]
        theta_col =                             m_theta[0]
    elif case == 2:
        theta_nc =                 m_theta[0] - m_theta[1]
        theta_col = (m_theta[0]+m_theta[1])/2 - theta_nc/2
    elif case == 3:
        theta_nc =             (m_theta[0] - m_theta[1])/2
        theta_col =            (m_theta[0] + m_theta[1])/2
    theta_col =                       180./np.pi*theta_col
    theta_nc =                    abs(180./np.pi*theta_nc)
    return(theta_col,theta_nc)

def setup_check(keyword,prefix,fname):

    #       CHECK FOR KEYWORD IN EITHER SETUP OR DIR, RETURN VAL

    pattern =   (keyword + "\s*=\s*['a-z]+|"
                + keyword + "\s*=\s*(\d+(.\d*)+(.\d*)E\+\d\d)|"
                + keyword + "\s*=\s*")

    with open(prefix + fname + '/setup.cfg','r') as fin:
        for num,line in enumerate(fin,1):
            s =     re.search(pattern,line)
            if s !=                   None:
                if s.group(2) !=      None:
                    val = float(s.group(1))
    return(val)

def spin_draw(ax,m_theta,pos1,pos2,colorin = "b"):
    
    #       DRAWS ARROWS ON THE GIVEN AXIS AX, ARROWS COORDINATES ARE SPECIFIED
    #       BY THE ANGLES M_THETA
    
    dx = np.sin(m_theta)
    dz = np.cos(m_theta)
    ax.arrow(pos1, 0, dx_col[1], dz_col[1], head_width=0.1, head_length=0.1,
                color = colorin)
    ax.arrow(pos2, 0, dx_col[0], dz_col[0], head_width=0.1, head_length=0.1,
                color = colorin)

def data_set_gen(inputname,dateset,key1 = "lambda",key2 = "scan",key3= "col"):

    #       GENERATES A STR SPINE TO INDICATE SOURCE DIR FOR DATA USED IN PLOT
    
    inputnamedummy = inputname
    newmatch = re.search(key1,inputname)
    newmatch2 = re.search(key2,inputname)
    newmatch3 = re.search(key3,inputname)
    string = ''
    if newmatch!= None:
        string = string + inputname[newmatch.start():newmatch.end()]
    if newmatch2!= None:
        string = string + "_" + inputname[newmatch2.start():]
    elif newmatch2 == None:
        while newmatch3!=None:
            newmatch3 = re.search('col',inputnamedummy)
            if newmatch3!=None:
                inputnamedummy = inputnamedummy[newmatch3.end():]
        string = string + '_col_' + inputnamedummy 
    retstr = "Data from: " + dateset + "/.." + string
    return(retstr)

def find_node(a,b,thresh):

    #

	mines = np.isclose(a,b,atol = thresh)
	return(mines)
weyl_node_idx = find_node(E_sq[2],E_sq[1],0.009)

def cont_plot(ax, zdata, diff = False, xdata = K[:2,:],  title ='Standardtitle',cbarlab = 'Energy in eV'
              , axlab = 'k', x_ticks = kx,y_ticks = ky, kx_lab = kx_lab, ky_lab = ky_lab):

    #PLOTS THE ENERGYLANDSCAPE ()

	sq_shape2= 1600
	sq_shape1= 3200
	E_sq = zdata.reshape(sq_shape1,sq_shape2)
	[X,Y] = xdata[:2,:].reshape(2,sq_shape1,sq_shape2)
	ax.set_aspect('equal')
	ax.set_title(title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="3%", pad=0.05)    
	norm = colors.Normalize(vmin = min(zdata), vmax = max(zdata))
	im3 = ax.contourf(X,Y,E_sq, cmap = cm.get_cmap(colmap),norm = norm
                      , levels = np.linspace(min(zdata),max(zdata),10))
	cbarticks = [min(zdata),np.mean(zdata),max(zdata)]
	if diff == True:
		cbarticks.append(0.0)
	cbarticklabels = ['{:1.1f}'.format(x) for x in cbarticks]
	cbar3 = fig.colorbar(im3,cax = cax,orientation = 'vertical',format = '%1.1f'
                        ,ticks = cbarticks)
	cbar3.set_label(cbarlab)

    #       GLOBAL VARS

findbool = False
files = 51
j = 0
peakexp = 3
colmap = 'jet'
datesets =       ["20190423","20190424","20190425","20190517","20190606","20190611"
                                                 ,"20190612","20190614","20190618"]
key =                                                               'anticol_theta'
rotation =                                                                'scancol'
rotation2 =                                                           'scancolflip'
seedfname =                                                   "path_rel_G-K-Kprime"
save_prefix =                          "/Users/kipp/STB/weyl_figures/Weyl_figures/"
coltypes =                                                         ['col','noncol']
func_col_types =                                                [1,0,2,3,2,2,2,2,2]
hc_syms =                                                         ['sym','antisym']
K_label =                               ["","$\Gamma$", "$K^{'}$", "K", "$\Gamma$"]
symdict =                                              {'sym':0,'antisym':files//2}
symstringdict =                      {'sym': 'symmetric','antisym':'antisymmetric'}
colthetalabeldict = {'col':'$\Theta_{col}$ [deg.]','noncol':'$\Theta_{nc}$ [deg.]'}
symhclabeldict =                       {'sym':r'$\sigma_{sym,xy}$ [$\frac{e^2}{h}$]'
                                  ,'antisym':r'$\sigma_{asym,xy}$ [$\frac{e^2}{h}$]'
                                        ,'bare':r"$\sigma_{xy}$ [$\frac{e^2}{h}$]"}
labeldict =        {'pwd':r"$\frac{\Delta\sigma}{\Delta\Theta}$",'peak':"$\sigma$"}
factordict =                           {'m_theta':180./np.pi,'lambda':1.,'t_so':1.}
prefixsold =         ["/Users/kipp/STB/output_jureca/" + d + "/" for d in datesets]
prefixs =                         ["/Data/ias-1/kipp/" + d + "/" for d in datesets]
save_strs =                                   ["datafrom_" + ds for ds in datesets]
inputname =                                  seedfname + '_' + key + '_' + rotation
inputname2 =                                seedfname + '_' + key + '_' + rotation2
coltypesdict =                             dict(list(zip(datesets,func_col_types)))
d = date.isoformat(date.today())
v = 0
norm = colors.Normalize(vmin = 0, vmax = 20)
lims = [-3.2,-2.8]
zoomlims = [-3.05,-2.90]
m,n,o = 1,6,18
prefix2 = prefixs[1]
prefix = prefixs[8]
dateset = datesets[8]
j = 2
key = "_mix"
inputname_null = "path_rel_G-K-Kprime_anticol_theta_col_{:02}".format(j) + key
inputname = "path_rel_G-K-Kprime_anticol_theta_scancol_{:02}".format(j) + key
inputname2 = "path_rel_G-K-Kprime_anticol_theta_scancolflip_{:02}".format(j) + key

    #       IMPORT

E_sq = np.load(prefix + inputname + '/band_E.npy')
E_null_sq = np.load(prefix + inputname_null + '/band_E.npy')
E_opp_sq = np.load(prefix + inputname2 + '/band_E.npy')

hc = np.load(prefix + inputname + '/hall_cond.npy')
hc_opp = np.load(prefix + inputname2 + '/hall_cond.npy')
hc_null = np.load(prefix + inputname_null + '/hall_cond.npy')
hc_E = np.load(prefix + inputname + '/hall_cond_E.npy')

m_theta = np.load(prefix + inputname + '/m_theta.npy')
m_theta_opp = np.load(prefix + inputname2 + '/m_theta.npy')
m_theta_null = np.load(prefix + inputname_null + '/m_theta.npy')

m_phi = np.load(prefix + inputname + '/m_phi.npy')
m_phi_opp = np.load(prefix + inputname2 + '/m_phi.npy')
m_phi_null = np.load(prefix + inputname_null + '/m_phi.npy')

lmda = setup_check('lambda',prefix,inputname + '/')
lmda_null = setup_check('lambda',prefix,inputname_null + '/')
t_so = setup_check('t_so',prefix,inputname + '/')
t_so_null = setup_check('t_so',prefix,inputname_null + '/')
temperature = setup_check('temperature',prefix,inputname + '/')

E = reshape_func(3200,1600,E_sq,4)
E_null = reshape_func(3200,1600,E_null_sq,4)
E_opp = reshape_func(3200,1600,E_opp_sq,4)
thetastr = "$\Theta_{cant}$" + " = {:1.2f} deg.".format(theta_nc) + "\n$\Theta_{col}$" + " = {:1.2f} deg.".format(theta_col)
lambdastr = "$\lambda_{xc}$" + " = {0:1.2f}".format(lmda)
t_sostr = "$t_{so}$" + " = {0:1.2f}".format(t_so)
temperaturestr = "$temperature$" + " = {0:1.2f} K".format(temperature)
textstr = thetastr + "\n" + temperaturestr + "\n" + lambdastr + "\n" + t_sostr

theta_col,theta_nc = col_type_func(m_theta,coltypesdict[dateset])
phi_col,phi_nc = col_type_func(m_phi,coltypesdict[dateset])
dx = np.sin(m_theta)
dz = np.cos(m_theta)
dx_opp = np.sin(m_theta_opp)
dz_opp = np.cos(m_theta_opp)
dx_col = np.sin(m_theta_null)
dz_col = np.cos(m_theta_null)
npoints_E = E.shape[1]
K = np.arange(npoints_E)

    #       PLOTTING

plt.style.use("/Users/kipp/Plotting/custom_style.mplstyle")

def band_compare_zoom(E,K,hc_fin,inputname,lims = lims,zoomlims = zoomlims):    #hc_fin = hc_fin_func(hc,hc_opp,"antisym")

    fig,ax = plt.subplots(nrows = 2, ncols = 3,figsize = (3*6.4,2*4.8),constrained_layout = True,sharex = "column",sharey = 'row')

    #       MAKE FAKEPLOTS

    for axes in [ax[0,0],ax[1,0]]:
        axes.plot([],[],color =cm.get_cmap(colmap)(norm(m)),label = "$E(+\Theta_{nc})$")
        axes.plot([],[],color =cm.get_cmap(colmap)(norm(n)),label = "$E(0)$")
        axes.plot([],[],color =cm.get_cmap(colmap)(norm(o)),label = "$E(-\Theta_{nc})$")
    for j in range(4):
        for axes in [ax[0,0],ax[1,0]]:
            axes.plot(K,E[j],color =cm.get_cmap(colmap)(norm(m)))
            axes.plot(K,E_opp[j],color =cm.get_cmap(colmap)(norm(n)))
            axes.plot(K,E_null[j],color =cm.get_cmap(colmap)(norm(o)))
    for axes in [ax[0,1],ax[1,1]]:
        axes.plot(hc_fin,hc_E,color =cm.get_cmap(colmap)(norm(m)))

    ax[0,2].set_axis_off()
    ylims = ax[0,2].get_ylim()
    xlims = ax[0,2].get_xlim()
    axins = inset_axes(ax[0,2], width=2.6, height=1.8,loc = 'upper right')
    axins.tick_params(tick1On = False,labelleft=False, labelbottom=False)
    axins.set_xlim(-1.8,2.8)
    axins.set_ylim(-1.2,1.2)
    color = cm.get_cmap(colmap)(norm(n))
    spin_draw(axins,m_theta,1.65,-0.65,color)
    spin_draw(axins,m_theta_col,1.65,-0.65,color)
    spin_draw(axins,m_theta_opp,1.65,-0.65,color)
    axes.text(xlims[0], ylims[0], textstr,size = 18,
                horizontalalignment='left',
                verticalalignment='bottom')

    #       CUSTOMIZING

    for axes in [ax[0,1],ax[1,1]]:
        axes.axhspan(ymin = zoomlims[0], ymax = zoomlims[1],facecolor = "b",alpha = 0.1)
        axes.set_xlabel(symhclabeldict["antisym"],fontsize = 18)
        axes.xaxis.set_major_locator(plt.MaxNLocator(3))
        axes.tick_params(labelsize = 18)
    for axes in [ax[0,0],ax[1,0]]:
        axes.axhspan(ymin = zoomlims[0], ymax = zoomlims[1],facecolor = "b",alpha = 0.1)
        axes.xaxis.set_major_locator(plt.MultipleLocator(npoints_E//3))
        axes.set_xticklabels(K_label)
        axes.legend(fontsize = 18)
        axes.set_ylabel("Energy [eV]",fontsize = 18)
        axes.margins(x=0, y=-0.0)
    ax[1,0].set_ylim(zoomlims)
    ax[0,1].set_ylim(lims)
    ax[1,2].set_axis_off()
    ylims = ax[1,2].get_ylim()
    xlims = ax[1,2].get_xlim()
    textstr = data_set_gen(inputname = inputname,dateset = dateset)
    ax[1,2].text(xlims[0], ylims[0],string,size = 18,
            horizontalalignment='left',
            verticalalignment='bottom',wrap = 'True')

    #       SAVING
    with PdfPages(save_prefix + d + '_' + inputname + '_' + 'bandscompare.pdf') as pdf:
        pdf.savefig(fig)
###################################################################################################################################################
colmap = 'Blues'
norm = colors.Normalize(vmin = -4, vmax = 4)
prefix = prefixs[8]
dateset = datesets[8]
filename = "path_rel_G-K-Kprime_anticol_theta_col_{:02}_large".format(2)
save_prefix = '/Users/kipp/STB/weyl_figures/Weyl_figures/'
save_fname = 'Bands_contour'
v = 0

s = date.isoformat(date.today())

E = np.load(prefix + filename + "/band_E.npy")
K = np.load(prefix + filename + "/band_k.npy")
rlat = np.load(prefix + filename + "/rez_lattice.npy")
E_shape = E.shape
K_shape = K.shape
sq_shape2 = 1600
sq_shape1 = 3200
E_sq = E.reshape(E_shape[0],sq_shape1,sq_shape2)
K_sq = K.reshape(K_shape[0],sq_shape1,sq_shape2)
E_line = reshape_func(3200,1600,E,4)
K_line = np.arange(E_line.shape[1])
m_theta = np.load(prefix + filename + '/m_theta.npy')


def contour_band_plotter(E,rlat,E_line,m_theta,rangex = 1,rangey = 2):
    theta_col,theta_nc = col_type_func(m_theta,coltypesdict[datesets[2]])
    dx = np.sin(m_theta)
    dz = np.cos(m_theta)
    kx = np.linspace(min(K[0,:]),max(K[0,:]),5)
    ky = np.linspace(min(K[1,:]),max(K[1,:]),5)
    kx_lab = ['{:1.1f}'.format(x) for x in kx]
    ky_lab = ['{:1.1f}'.format(x) for x in ky]
    vx = rlat[:,0]
    vy = rlat[:,1]
    G_point = vx + vy
    M_point = 0.5*vx + 0.5*vy
    K_point = 1./3*vx + 2./3*vy
    K_prime_point = 2./3*vx + 1./3*vy
    path_verts = np.array([[0,0],K_prime_point,K_point])
    path = patches.Polygon(path_verts,fill = None, edgecolor = 'k',linewidth = 0.5, linestyle = ':')
    
    fig,axis = plt.subplots(nrows = 2, ncols = 3,figsize = (3*6.4/2,2*4.8),constrained_layout = True)
    ############plt.subplots_adjust(wspace = 0.8,hspace = 0.0)
    inputvec = 0.5/np.cos(np.pi/6)*np.dot(ucv.R(0.5),vy)
    hexcoords = ucv.hexcoords(inputvec)
    axes = [axis[0,0],axis[0,1],axis[1,0],axis[1,1]]
    for j in range(4):
	    cont_plot(ax = axes[j],zdata = E[j],title = 'Band_{}'.format(j+1))
	    path = patches.Polygon(path_verts + vx + vy,fill = None, edgecolor = 'r',lw = 1.0)
	    axes[j].add_patch(path)
	    for l in range(1,rangex+1):
		    for h in range(1,rangey+1):
			    hexagon = patches.Polygon(hexcoords +l*vx + h*vy ,fill = None, lw= 0.5, edgecolor = 'k')
			    axes[j].add_artist(hexagon)
    path = patches.Polygon(path_verts + vx + vy,fill = None, edgecolor = 'r',lw = 1.0)
    axis[0,2].add_patch(path)
    for l in range(1,rangex+1):
	    for h in range(1,rangey+1):
		    hexagon = patches.Polygon(hexcoords +l*vx + h*vy ,fill = None, lw= 0.5, edgecolor = 'k')
		    axis[0,2].add_artist(hexagon)
    cont_plot(ax = axis[0,2],diff = True, zdata = E[2]-E[1],title = "E_{0} - E_{1}".format(3,2))

    lmda = setup_check('lambda',prefix,filename + '/')
    t_so = setup_check('t_so',prefix,filename + '/')
    temperature = setup_check('temperature',prefix,filename+ '/')
    thetastr = "$\Theta_{cant}$" + " = {:1.2f} deg.".format(theta_nc) + "\n$\Theta_{col}$" + " = {:1.2f} deg.".format(theta_col)
    lambdastr = "$\lambda_{xc}$" + " = {0:1.2f}".format(lmda)
    t_sostr = "$t_{so}$" + " = {0:1.2f}".format(t_so)
    temperaturestr = "$temperature$" + " = {0:1.2f} K".format(temperature)
    textstr = thetastr + "\n" + temperaturestr + "\n" + lambdastr + "\n" + t_sostr
    datastr = data_set_gen(inputname = inputname,dateset = dateset)
    textstr = textstr + datastr

    ylims = axis[1,2].get_ylim()
    xlims = axis[1,2].get_xlim()
    axins = inset_axes(axis[1,2], width=2.0, height=1.2,loc = 'upper right')
    axins.tick_params(tick1On = False,labelleft=False, labelbottom=False)
    axins.set_xlim(-1.8,2.8)
    axins.set_ylim(-1.2,1.2)
    spin_draw(axins,m_theta,1.65,-0.65,colorin = "b")
    axis[1,2].text(xlims[0], ylims[0], textstr,size = 18,
             horizontalalignment='left',
             verticalalignment='bottom')
    axis[1,2].set_axis_off()
    axis[1,0].set_xlabel("$k_x$")
    axis[1,1].set_xlabel("$k_x$")
    axis[0,0].set_ylabel("$k_y$")
    axis[1,0].set_ylabel("$k_y$")
    with PdfPages(save_prefix + d + '_' + dateset + '_' + filename + '_' + save_fname + '.pdf') as pdf:
	    pdf.savefig(fig)