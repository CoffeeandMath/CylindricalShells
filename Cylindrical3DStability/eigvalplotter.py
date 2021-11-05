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


ax.plot(eigrenorm(numpy_array[0]))
ax.plot(eigrenorm(numpy_array[1]))
ax.plot(eigrenorm(numpy_array[2]))
ax.plot(eigrenorm(numpy_array[3]))
ax.plot(eigrenorm(numpy_array[4]))


plt.show()