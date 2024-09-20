import numpy as np
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
import matplotlib.colors as colors
import sys
import os
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
from numpy import inf

t = 0.95
a = np.loadtxt(f'data-{t:.4f}')
g1 = np.loadtxt(f'interface-{t:.2f}.dat')
x = a[:,0]
y = a[:,1]
f = a[:,2]
txx = a[:,9]
tyy = a[:,11]
tqq = a[:,12]
tp = np.log(txx+tyy+tqq)/np.log(10)
tp[np.isnan(tp)] = -10
tp[tp==-inf] = -10
tp[np.where(f==0)] = 100
solid = a[:,17]
nsolid = solid
nsolid[nsolid == 0.5] = 1.0
nsolid[nsolid == -1] = 0.0

N = 2048
xmin = -4
xmax = 4
ymin = 0
ymax = 8

X,Y = np.meshgrid(np.linspace(xmin,xmax,N,endpoint=False),np.linspace(ymin,ymax,N,endpoint=False))

Z = griddata((x, y), tp, (X, Y),method='linear')
zd1r = np.transpose(Z)

Z = griddata((x, y), nsolid, (X, Y),method='linear')
zd1rr = np.transpose(Z)

cmap = cm.get_cmap('hot')
norm = plt.Normalize(vmin=-1.3,vmax=0.1)
print(np.min(zd1r),np.max(zd1r),flush=True)
x_flat = X.flatten()
z_flat = Y.flatten()
magnitude_flat = np.transpose(zd1r).flatten()

colors = cmap(norm(magnitude_flat))

alpha = np.ones_like(colors[:,0])
is_white = np.all(colors[:,:3] >= 0.99,axis=1)
alpha[is_white] = 0 

output_data = np.column_stack((x_flat,z_flat,magnitude_flat,colors[:,:3],alpha))
output_file = 'blender_data.txt'
np.savetxt(output_file,output_data)

t = 0.95
a = np.loadtxt(f'data-{t:.4f}')
h = np.loadtxt(f'index-{t:.4f}')
x = a[:,0]
y = a[:,1]
f = a[:,2]
fi = h[:,6]
fi[f == 0] = 0


N = 2048
xmin = -4
xmax = 4
ymin = 0
ymax = 8

X,Y = np.meshgrid(np.linspace(xmin,xmax,N,endpoint=False),np.linspace(ymin,ymax,N,endpoint=False))

Z = griddata((x, y), fi, (X, Y),method='linear')
zd1r = np.transpose(Z)

Z = griddata((x, y), f, (X, Y),method='linear')
zdf = np.transpose(Z)

cmap = cm.get_cmap('RdBu')
norm = plt.Normalize(vmin=-1.,vmax=1.)

x_flat = X.flatten()
z_flat = Y.flatten()
magnitude_flat = np.transpose(zd1r).flatten()
f_flat = np.transpose(zdf).flatten()

colors = cmap(norm(magnitude_flat))

alpha = np.ones_like(colors[:,0])
#is_white = np.all(f_flat == 0.,axis=0)
alpha[f_flat < 0.5] = 0 

output_data = np.column_stack((x_flat,z_flat,magnitude_flat,colors[:,:3],alpha,f_flat))
output_file = 'blender_index_data.txt'
np.savetxt(output_file,output_data)
