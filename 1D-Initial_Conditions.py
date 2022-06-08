import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack as fft

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
u=np.zeros(K) #f(t=n)

if random==1:
    a=np.zeros(K)*1j
    
    for i in range(int(K/2)):
        ran=np.random.rand(1)*10**(-1)
        phi=np.random.rand(1)*2*np.pi
        a[i]=ran*np.exp(1j*phi)
        a[-i]=ran*np.exp(-1j*phi)
    
    a=fft.fftshift(a)
    
    u0=np.real(fft.ifft(a))

else:
        
    A=1
    sigma=0.4
    
    u0=A*np.exp(-(x-L/2)**2/sigma**2)
    a=fft.fft(u0)


plt.plot(k,abs(a),'o')
plt.xlabel(r"$k$")
plt.ylabel(r"$|a_k|$")
plt.tick_params(axis='both',direction='in',width=1,length=6)
plt.show()

plt.plot(u0)
plt.xlabel(r"$x$")
plt.ylabel(r"$u$")
plt.tick_params(axis='both',direction='in',width=1,length=6)
plt.show()

np.savetxt(r'u0.txt', u0)
