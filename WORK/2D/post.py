#!/usr/bin/python
import scipy as sp
import matplotlib.pyplot as plt
data1 = sp.genfromtxt('data/density_087.d')
x1 = data1[:,0]
y1 = data1[:-1,1]
plt.plot(x1,y1)
plt.savefig('data/density.png')
plt.close()

data2 = sp.genfromtxt('data/U_Velocity_087.d')
x2 = data2[:,0]
y2 = data2[:-1,1]
plt.plot(x2,y2)
plt.savefig('data/U_Velocity.png')
plt.close()

data3 = sp.genfromtxt('data/pressure_087.d')
x3 = data3[:,0]
y3 = data3[:-1,1]
plt.plot(x3,y3)
plt.savefig('data/pressure.png')
plt.close()

data4 = sp.genfromtxt('data/Mach_087.d')
x4 = data4[:,0]
y4 = data4[:-1,1]
plt.plot(x4,y4)
plt.savefig('data/Mach.png')
plt.close()

data = sp.genfromtxt('residual.d')
x = data[:,0]
y = data[:-1,1]
plt.yscale('log')
plt.plot(x,y)
plt.savefig('redidual.png')
plt.close()
