#!/usr/bin/python
import scipy as sp
import matplotlib.pyplot as plt

data = sp.genfromtxt('residual.dat',delimiter = '\t')
x = data[:,0]
y = data[:,1]
plt.plot(x,y)
plt.grid()
plt.show()
