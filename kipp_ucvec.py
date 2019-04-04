import matplotlib.pyplot as plt
from matplotlib import patches as pa
from matplotlib import cm,colors
import numpy as np

j = 1
prefix = '/Users/kipp/STB/output/weyl_theta_{}/'.format(j)
saveprefix = '/Users/kipp/STB/weyl_figures/'
lat = np.load(prefix + 'lattice.npy')
rlat = np.load(prefix + 'rez_lattice.npy')
path = np.load(prefix + 'band_k.npy')
vx = lat[:,0]
vy = lat[:,1]
kx = rlat[:,0]
ky = rlat[:,1]
R =lambda j: np.array([[np.cos(np.pi/3*j),np.sin(np.pi/3*j)],[-np.sin(np.pi/3*j),np.cos(np.pi/3*j)]])
vectors = [np.array([[0,0],x]) for x in [kx,ky]]
vec_lab = ['$k_1$','$k_2$']
M = 0.5*vectors[0]+0.5*vectors[1]
K = 1./3*vectors[0]+2./3*vectors[1]
Kp = 2./3*vectors[0] + 1./3*vectors[1]

def hexcoords(vector):	
	R =lambda j: np.array([[np.cos(np.pi/3*j),np.sin(np.pi/3*j)],[-np.sin(np.pi/3*j),np.cos(np.pi/3*j)]])
	hexcoords =np.array([np.dot(R(j),vector) for j in range(6)])
	return(hexcoords)

	#ky links unit cell centers, so ky/2. also we have to rotate ky by 30 degrees because the ky joins the centers crossing a side of the hexagon not an edge
'''
inputvector = 0.5/np.cos(np.pi/6)*np.dot(R(0.5),ky)
hexcoords = hexcoords(inputvector)
hexagon = pa.Polygon(hexcoords,fill=None,edgecolor = 'k',linewidth = 0.5)
hexagon2 = pa.Polygon(hexcoords + ky,fill= None,edgecolor = 'r',linewidth = 0.5)

norm = colors.Normalize(vmin = 0,vmax = 10)
fig,ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
ax.add_patch(hexagon)
ax.add_patch(hexagon2)
for j,x in enumerate(vectors):
	plt.plot(x[:,0],x[:,1],color = cm.get_cmap('hsv')(norm(j)),label =vec_lab[j])
annolabel = ['$K$',"$K^{'}$"]
for j,x in enumerate(hexcoords):
	ax.annotate(annolabel[j%2],xy = x,xycoords = 'data')
ax.annotate('$\Gamma$',xy = [0,0],xycoords = 'data',
				verticalalignment = 'bottom',horizontalalignment = 'left')
plt.plot(path[0],path[1],label = '$\Gamma$MK$\Gamma$')
plt.plot(M[:,0],M[:,1],'ko')
plt.plot(K[:,0],K[:,1],'ko')
plt.plot(Kp[:,0],Kp[:,1],'ko')
plt.legend()
plt.savefig(saveprefix + 'unitcell_vectors',format = 'pdf')
#plt.show()
'''
