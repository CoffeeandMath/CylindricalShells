import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numpy as np
from numpy import random 
import pandas as pd


index = 10

vrdata = pd.read_csv('Cylindrical3DStability/build/unstablemodes/vr_' + str(index) + '.csv')
vrsorted = vrdata.sort_values('S')
vthetadata = pd.read_csv('Cylindrical3DStability/build/unstablemodes/vtheta_'  + str(index) + '.csv')
vthetasorted = vthetadata.sort_values('S')
vzdata = pd.read_csv('Cylindrical3DStability/build/unstablemodes/vz_' + str(index) + '.csv')
vzsorted = vzdata.sort_values('S')

vd = {'S':vrsorted['S'],'r':vrsorted['val'],'theta':vthetasorted['val'],'z':vzsorted['val']}
vn = pd.DataFrame(data=vd)

wrdata = pd.read_csv('Cylindrical3DStability/build/unstablemodes/wr_' + str(index) + '.csv')
wrsorted = wrdata.sort_values('S')
wthetadata = pd.read_csv('Cylindrical3DStability/build/unstablemodes/wtheta_'  + str(index) + '.csv')
wthetasorted = wthetadata.sort_values('S')
wzdata = pd.read_csv('Cylindrical3DStability/build/unstablemodes/wz_' + str(index) + '.csv')
wzsorted = wzdata.sort_values('S')

wd = {'S':wrsorted['S'],'r':wrsorted['val'],'theta':wthetasorted['val'],'z':wzsorted['val']}
wn = pd.DataFrame(data=wd)

ntheta = 300
xrdf = pd.read_csv('CylindricalSystem/build/solutions/r_values_' + str(index) + '.csv')
xrdfsorted = xrdf.sort_values('S_values')
xzdf = pd.read_csv('CylindricalSystem/build/solutions/z_values_' + str(index) + '.csv')
xzdfsorted = xzdf.sort_values('S_values')
th = np.linspace(0, np.pi*2, ntheta)
xr = xrdfsorted['r_values']
xz = xzdfsorted['z_values']


nf = 1.;
uxn = np.outer(vn['r'],np.multiply(np.cos(th),np.cos(nf*th))) - np.outer(vn['theta'],np.multiply(np.sin(th),np.cos(nf*th))) + np.outer(wn['r'],np.multiply(np.cos(th),np.sin(nf*th))) - np.outer(wn['theta'],np.multiply(np.sin(th),np.sin(nf*th)))
uyn = np.outer(vn['r'],np.multiply(np.sin(th),np.cos(nf*th))) + np.outer(vn['theta'],np.multiply(np.cos(th),np.cos(nf*th))) + np.outer(wn['r'],np.multiply(np.sin(th),np.sin(nf*th))) + np.outer(wn['theta'],np.multiply(np.cos(th),np.sin(nf*th)))
uzn = np.outer(vn['z'],np.cos(nf*th)) + np.outer(wn['z'],np.sin(nf*th))

defmag = 2000.0
xxn = np.outer(xr, np.cos(th)) + defmag*uxn
xyn = np.outer(xr, np.sin(th)) + defmag*uyn
xzn = np.outer(xz,np.linspace(1.,1.,th.size)) + defmag*uzn

x0xn = np.outer(xr, np.cos(th)) 
x0yn = np.outer(xr, np.sin(th)) 
x0zn = np.outer(xz,np.linspace(1.,1.,th.size))
ax2 = [None]*2
fig = plt.figure(figsize=plt.figaspect(0.5))
ax2[0] = fig.add_subplot(121)
ax2[1] = fig.add_subplot(122,projection='3d')
ax2[0].plot(vn['S'],vn['r'],label='r')
ax2[0].plot(vn['S'],vn['theta'],label = 'theta')
ax2[0].plot(vn['S'],vn['z'], label = 'z')
ax2[0].legend()

ax2[1].plot_surface(xxn,xyn,xzn,linewidth=10, antialiased=True)
ax2[1].plot_wireframe(x0xn,x0yn,x0zn+1.)
axlim = 1.5
ax2[1].set_xlim3d([-axlim,axlim])
ax2[1].set_ylim3d([-axlim,axlim])
ax2[1].set_zlim3d([-axlim,axlim])
plt.show()