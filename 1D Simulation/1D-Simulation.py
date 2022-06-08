import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
import pandas as pd
from scipy import fftpack as fft
import sys

equation = 'WE'
# 'WE' = Wave Equation
# 'KG' = Klein-Gordon Equation
# 'SG' = Sine-Gordon Equation
# 'V4' = Quartic Interaction Equation
# 'uux' = u ∂u term
# 'uux2' = u (∂u)² term
# 'u2ux2' = u² (∂u)² term


ntau=5 #duration of the simulation in units of 𝜏=c·L
K=2**6 #Number of space steps and frequencies

gif = 1
# 1 to turn on gif
# 0 to turn off gif

# if gif is different than 0 or 1 gif will be turn off. 

# if gif turned is chosen the duration of the simulation would be limited to t=30𝜏

#Constants
c=1
L=2*np.pi

m=0.1
lam=0.5
alpha=0.01
beta=0.001
gamma=0.001

#Discretitzation
dx=L/K
krange=np.arange(-K/2,K/2)
k=fft.fftshift(krange)

tau=c*L
if gif == 1 and ntau>30: ntau=30

t0=0.
t_bound=tau*ntau

d=5
dt=dx/c/d

N=int((t_bound-t0)/dt)
print('N =',N)

#Set up variables
x=np.linspace(0,L-dx,K)
u=np.zeros(K) #f(t=n)
U=np.zeros(K) #f(t=n+1)
v=np.zeros(K)
V=np.zeros(K)
w=np.zeros(K)
W=np.zeros(K)

t=np.linspace(t0,t_bound,N)


#Define Functions

#Derivative
def D(u):
    a=fft.fft(u)
    b=k*a*1j    
    ux=np.real(fft.ifft(b))
    return ux

#Differential Equation
def f(a1,b1,c1):
    if equation=='WE': nonlinearterm=0
    if equation=='KG': nonlinearterm=-m**2*a1
    if equation=='SG': nonlinearterm=-m**2*np.sin(np.pi/2*a1)
    if equation=='V4': nonlinearterm=-m**2*a1-4*lam*a1**3
    if equation=='uux': nonlinearterm=-alpha*a1*c1
    if equation=='uux2': nonlinearterm=-beta*a1*c1**2
    if equation=='u2ux2': nonlinearterm=-gamma*a1**2*c1**2
 
    a2=b1
    b2=c**2*D(c1)+nonlinearterm
    c2=D(b1)
    return a2,b2,c2

#Runge-Kutta 4
def RK4(a1,b1,c1):
    K1=f(a1,b1,c1)
    K2=f(a1+dt/2*K1[0],b1+dt/2*K1[1],c1+dt/2*K1[2])
    K3=f(a1+dt/2*K2[0],b1+dt/2*K2[1],c1+dt/2*K2[2])
    K4=f(a1+dt*K3[0],b1+dt*K3[1],c1+dt*K3[2])
    return a1+dt/6*(K1[0]+2*K2[0]+2*K3[0]+K4[0]),b1+dt/6*(K1[1]+2*K2[1]+2*K3[1]+K4[1]),c1+dt/6*(K1[2]+2*K2[2]+2*K3[2]+K4[2])

#Initial Conditions
u0 = np.array(pd.read_csv('u0.txt',header=None)[0])
u=u0
w=D(u)

#Set up gif
if gif==1:
    fig, ax = plt.subplots()
    camera = Camera(fig)

#Check if equation variable is well defined
if equation!='WE' and equation!='KG' and equation!='SG' and equation!='V4' and equation!='uux' and equation!='uux2' and equation!='u2ux2':
    sys.exit('Error: equation value is not valid')

#Define residue to make 60 frames
r=N//59

#Time loop
for n in range(N):
    if n%r==0:
        print("{}/{}".format(n,N))    
        
        #Plot gif frame
        if gif==1:
            plt.plot(x,u,label="t={}".format(n),color='C0')
            plt.xlim(0,L)
        
            camera.snap()

    #Calculate variables of the next time-step
    U,V,W=RK4(u,v,w)
    
    #Break the loop if the simulation diverge
    if np.isnan(U[0])==True:
        break

    #Replace old variables with new variables
    u=U
    v=V
    w=W

    
print("Finished")

#Process gif
if gif==1:
    animation = camera.animate()
    animation.save('1D_Simulation.gif')
    plt.show()

#Plot Power Spectrum
a=fft.fft(u)

plt.plot(k[:int(len(k)/2)],abs(a[:int(len(k)/2)]),'--')
plt.xlabel(r"$k$")
plt.ylabel(r"$|a_k|$")
plt.yscale('log')
plt.tick_params(axis='both',direction='in',width=1,length=6)
plt.show()

#Save Power Spectrum
np.savetxt(r'PS_k.txt', k, fmt='%d')
np.savetxt(r'PS_'+equation+'.txt', abs(a))
