import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numpy as np
from numpy import random 
import pandas as pd


index = 666

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
th = np.linspace(0, np.pi*2, ntheta)
xr = np.linspace(0.1,1.,vn['r'].size)
xz = np.linspace(0.,0.,vn['r'].size)

uxn = np.outer(vn['r'],np.multiply(np.cos(th),np.cos(th))) - np.outer(vn['theta'],np.multiply(np.sin(th),np.cos(th))) + np.outer(wn['r'],np.multiply(np.cos(th),np.sin(th))) - np.outer(wn['theta'],np.multiply(np.sin(th),np.sin(th)))
uyn = np.outer(vn['r'],np.multiply(np.sin(th),np.cos(th))) + np.outer(vn['theta'],np.multiply(np.cos(th),np.cos(th))) + np.outer(wn['r'],np.multiply(np.sin(th),np.sin(th))) + np.outer(wn['theta'],np.multiply(np.cos(th),np.sin(th)))
uzn = np.outer(vn['z'],np.cos(th)) + np.outer(wn['z'],np.sin(th))

defmag = 2.0
xxn = np.outer(xr, np.cos(th)) + defmag*uxn
xyn = np.outer(xr, np.sin(th)) + defmag*uyn
xzn = np.outer(xz,np.linspace(1.,1.,th.size)) + defmag*uzn

fig2, ax2 = plt.subplots(1,2)
ax2[1] = fig2.add_subplot(122,projection='3d')
ax2[0].plot(vn['S'],vn['r'],label='r')
ax2[0].plot(vn['S'],vn['theta'],label = 'theta')
ax2[0].plot(vn['S'],vn['z'], label = 'z')
ax2[0].legend()

ax2[1].plot_surface(xxn,xyn,xzn,linewidth=10, antialiased=True)
ax2[1].set_xlim3d([-2,2])
ax2[1].set_ylim3d([-2,2])
ax2[1].set_zlim3d([-2,2])
plt.show()