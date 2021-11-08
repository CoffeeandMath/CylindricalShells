import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

plt.style.use('default')
fig, ax = plt.subplots(nrows=1,ncols=1)

evalsf = open("build/eigenvalues.csv")
numpy_array = np.loadtxt(evalsf,delimiter=",")

normalized = True

def eigrenorm(array):
	if (normalized):
		return array/array[0]
	else:
		return array


ax.plot(eigrenorm(numpy_array[0]),label='n = 1')
ax.plot(eigrenorm(numpy_array[1]),label='n = 2')
ax.plot(eigrenorm(numpy_array[2]),label='n = 3')
ax.plot(eigrenorm(numpy_array[3]),label='n = 4')
ax.plot(eigrenorm(numpy_array[4]),label='n = 5')
ax.legend(loc='upper right')

plt.show()