import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack as fft
import pandas as pd

random = 0
# 1 to turn on weak random initial conditions
# 0 to take a Gaussian as initial condition
# if random = 0 or is different than 1 it will create a Gaussian

c=1
L=2*np.pi

K=2**6

dx=L/K
krange=np.arange(-K/2,K/2)
k=fft.fftshift(krange)

x=np.linspace(0,L-dx,K)
y=np.linspace(0,L-dx,K)
u=np.zeros(K) #f(t=n)

if random==1:
    
    a=np.zeros((K,K))*1j
    for i in range(int(K)):
        for j in range(int(K)):
            ran=np.random.rand(1)*10**(-1)
            phi=np.random.rand(1)*2*np.pi
            a[i,j]=ran[0]*np.exp(1j*phi[0])
            
    u0=np.real(fft.ifft(a))

else:
    A=1
    sigma=0.7
    xx,yy = np.meshgrid(x,y)
    u0=A*np.e**(-(xx-L/2)**2/sigma**2)*np.e**(-(yy-L/2)**2/sigma**2)

a=fft.fft2(u0)
plt.pcolormesh(krange,krange,abs(a),cmap='inferno')
plt.colorbar()
plt.show()
plt.pcolormesh(x,y,u0,cmap='inferno')
plt.colorbar()
plt.show()

with open('u0.txt','w') as f:
    for i in range(K):
        for j in range(K):
            f.write("{}".format(u0[i,j]))
            if j!=K-1: f.write(" ")
        if i!=K-1: f.write("\n")
           
