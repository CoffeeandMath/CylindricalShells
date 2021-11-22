import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from numpy import random 
import pandas as pd

evalsymdata = pd.read_csv('CylindricalSystem/build/evals.csv').values
num_rows, num_cols = evalsymdata.shape
color = 'black'
normalized = 1


fig, ax1 = plt.subplots()

def nanv(v):
	vout = [None]*len(v)
	for i in range(len(v)):
		if np.fabs(v[i]) > 1.e-10:
			vout[i] = v[i]
		else:
			vout[i] = np.nan
	return vout



def normalize(v,tf):
	if tf==1:
		vtemp = [None]*len(v)
		for i in range(len(v)):
			vtemp[i] = v[i]/v[0]

		return vtemp 
	vtemp = v.copy()
	return vtemp


def getcolumn(Mat,k):
	earray = [None]*len(Mat);
	for i in range(len(Mat)):
		earray[i] = Mat[i][k]

	return earray



for i in range(0,num_cols):
	ax1.plot(normalize(evalsymdata,normalized), color = color)


ax1.axhline(y=0.0,color = 'blue',linestyle='--')
#ax1.set_ylim([-1.0e-6,1.0e-6])
#Plotting the asymmetric eigenvalues
normalized = 0
nevals = 5
fig2, ax2 = plt.subplots(1,nevals)
evalasymdata = [None] * nevals;



for i in range(1,nevals+1):
	evtemp = pd.read_csv('Cylindrical3DStability/build/evals'+str(i)+'.csv').values
	jstart = 1
	ax2[i-1].title.set_text(str(i))
	if i==1:
		jstart = 1

	for j in range(jstart,len(evtemp[0])+1):
		coldata = getcolumn(evtemp,j-1)
		ax2[i-1].plot(nanv(normalize(coldata,normalized)), color = color)

for i in range(nevals):
	ax2[i].axhline(y=0.0,color = 'blue',linestyle='--')

for i in range(nevals):
	ax2[i].set_xlim([0.,num_rows])
ax2[0].set_ylim([-1.0e-3,1.0e-3])
ax2[1].set_ylim([-0.005,0.02])
ax2[2].set_ylim([-0.05,0.25])
ax2[3].set_ylim([-0.05,0.4])
ax2[4].set_ylim([-0.05 , 0.5])
fig2.suptitle("Fourier Mode Stability")
plt.show()