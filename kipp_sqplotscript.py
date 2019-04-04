import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,ticker,colors,patches
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages
import kipp_ucvec as ucv
from datetime import date
prefix = '/Users/kipp/STB/output/sq_fill_anticol_anticol_theta_00/'
prefixline = '/Users/kipp/STB/output/path_rel_G-K-Kprime_anticol_anticol_theta_00/'
filenameE = 'band_E.npy'
filenameEline = 'band_E.npy'
filenameKline = 'band_k.npy'
filenameK = 'band_k.npy'
filenamelat = 'rez_lattice.npy'
save_prefix = '/Users/kipp/STB/weyl_figures/Weyl_figures/'
save_fname = 'Bands_contour_anticol'
d = '2019-02-27'
v = 0

'''today = date.isoformat(date.today())
today_old  = np.load(prefix + 'today.npy')
print(today_old,today)
if today != today_old:
        np.save(prefix + 'today',today)
if today==today_old:
        version = np.load(prefix + 'version.npy') + 1
        print(version)
        np.save(prefix + 'version.npy',version)

prefix = '/Users/kipp/STB/output/'
filenameE = 'sqfill_1uc_thetapihalf/band_E.npy'
filenameEline = 'weyl_theta_sqline/band_E.npy'
filenameKline = 'weyl_theta_sqline/band_k.npy'
filenameK = 'sqfill_1uc_thetapihalf/band_k.npy'
filenamelat = 'sqfill_1uc_thetapihalf/rez_lattice.npy'
save_prefix = '/Users/kipp/STB/weyl_figures/Weyl_figures/'
save_fname = 'Bands_contour'
'''
E = np.load(prefix + filenameE)
Eline = np.load(prefixline + filenameEline)
Kline = np.load(prefixline + filenameKline)
K = np.load(prefix + filenameK)
rlat = np.load(prefix + filenamelat)

E_shape = E.shape
K_shape = K.shape
sq_shape = int(np.sqrt(E_shape[1]))

E_sq = E.reshape(E_shape[0],sq_shape,sq_shape)

k_ticks = np.linspace(0,sq_shape,5)
kx = np.linspace(min(K[0,:]),max(K[0,:]),5)
ky = np.linspace(min(K[1,:]),max(K[1,:]),5)
kx_lab = ['{:1.1f}'.format(x) for x in kx]
ky_lab = ['{:1.1f}'.format(x) for x in ky]

vx = rlat[:,0]
vy = rlat[:,1]
M_point = 0.5*vx + 0.5*vy
K_point = 1./3*vx + 2./3*vy
K_prime_point = 2./3*vx + 1./3*vy
path_verts = np.array([[0,0],K_prime_point,K_point])
path = patches.Polygon(path_verts,fill = None, edgecolor = 'k',linewidth = 0.5, linestyle = ':')

def find_node(a,b,thresh):
	mines = np.isclose(a,b,atol = thresh)
	return(mines)
weyl_node_idx = find_node(E_sq[2],E_sq[1],0.009)

def cont_plot(ax, zdata, xdata = K[:2,:],  title ='Standardtitle',cbarlab = 'Energy in eV', axlab = 'k', x_ticks = kx,y_ticks = ky, kx_lab = kx_lab, ky_lab = ky_lab):
	sq_shape = int(np.sqrt(E_shape[1]))
	E_sq = zdata.reshape(sq_shape,sq_shape)
	[X,Y] = xdata[:2,:].reshape(2,sq_shape,sq_shape)
	ax.set_aspect('equal')
	xaxlab = '$' + axlab + '_x$'
	yaxlab = '$' + axlab + '_y$'
	
	ax.set_xlabel(xaxlab)
	ax.set_ylabel(yaxlab)
	ax.set_xticks(x_ticks)
	ax.set_yticks(y_ticks)
	ax.set_title(title)
	im3 = ax.contourf(X,Y,E_sq, cmap = cm.get_cmap('RdBu'),levels = 10)
	cbarticks = [min(zdata),np.mean(zdata),max(zdata)]
	cbarticklabels = ['{:1.1f}'.format(x) for x in cbarticks]
	cbar3 = fig.colorbar(im3,ax = ax,orientation = 'horizontal',format = '%1.1f',ticks = cbarticks)
	cbar3.set_label(cbarlab)

fig,axis = plt.subplots(nrows = 2, ncols = 4,figsize = (12,12),constrained_layout = True)
inputvec = 0.5/np.cos(np.pi/6)*np.dot(ucv.R(0.5),vy)
hexcoords = ucv.hexcoords(inputvec)

for j in range(4):
	print(j)
	if j>0:
		axis[1,j].remove()
	cont_plot(ax = axis[0,j],zdata = E[j],title = 'Band_{}'.format(j+1))
	path = patches.Polygon(path_verts + vx + vy,fill = None, edgecolor = 'g',lw = 1.0)
	axis[0,j].add_patch(path)
	for l in range(1,3):
		for h in range(1,3):
			hexagon = patches.Polygon(hexcoords +l*vx + h*vy ,fill = None, lw= 0.5, edgecolor = 'k')
			axis[0,j].add_artist(hexagon)
			
cont_plot(ax = axis[1,0],zdata = E[2]-E[1],title = "E_{0} - E_{1}".format(3,2))

'''for l in range(1,3):
	for h in range(1,3):
		hexagon = patches.Polygon(hexcoords +l*vx + h*vy ,fill = None, lw= 0.5, edgecolor = 'k')
		axis[1,0].add_artist(hexagon)
'''
#path = patches.Polygon(path_verts + vx + vy,fill = None, edgecolor = 'g',lw = 1.0)

#[x_circ,y_circ] =  K[:2,:].reshape(2,sq_shape,sq_shape)
#xy = np.array([x_circ[weyl_node_idx],y_circ[weyl_node_idx]]).transpose()
#circles = [patches.Ellipse(B,width = 1.,height = 1., edgecolor = 'g', lw = 1., fill = None) for B in xy]
#axis[1,0].add_patch(path)
#for c in circles:
#	axis[1,0].add_artist(c)
#cont_plot(ax = axis[1,1],zdata = E[1]-E[0],title = "E_{0} - E_{1}".format(2,1))
'''
norm = colors.Normalize(vmin = 0,vmax = 4)
K_plot = np.linspace(0,max(Kline[0]),Kline.shape[1])
#print(K_.shape,Eline.shape)
for j in range(Eline.shape[0]):
	axis[1,0].plot(K_plot,Eline[j],color = cm.get_cmap('plasma')(norm(j)))
axis[1,0].set_title("Bandstructure along $\Gamma$,$K^{'}$,K")
axis[1,0].set_xticks(max(Kline[0])/3*np.array([0,1,2,3]))
axis[1,0].set_xticklabels(['$\Gamma$',"$K^{'}$",'K','$\Gamma$'])
'''
plt.show()
#with PdfPages(save_prefix + d + '_' + str(v) + '_'  + save_fname + '_anticol_theta_0.pdf') as pdf:
#	pdf.savefig(fig)	
#plt.show()

