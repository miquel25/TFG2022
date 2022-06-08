import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
from scipy import fftpack as fft
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import copy
from scipy import stats
import sys

equation = 'WE'
# 'WE' = Wave Equation
# 'KG' = Klein-Gordon Equation
# 'SG' = Sine-Gordon Equation
# 'V4' = Quartic Interaction Equation
# 'uux' = u ∇u term
# 'uux2' = u (∇u)² term
# 'u2ux2' = u² (∇u)² term

# if equation is different than any of these options it will simulate the Wave Equation

ntau=3 #duration of the simulation in units of 𝜏=c·L

gifDimensions = 3
# 2 to turn on gif in 2D
# 3 to turn off gif in 3D

# if gif is different than 2 or 3 gif will be turn off

# if gif turned is chosen the duration of the simulation would be limited to t=30𝜏

#Constants
c=1
L=2*np.pi

m=0.1
lam=0.001
alpha=0.5
beta=0.5
gamma=0.5

#Discretitzation
K=2**6

dx=L/K
dy=dx
krange=np.arange(-K/2,K/2)
kx=fft.fftshift(krange)
ky=fft.fftshift(krange)
k2d=np.meshgrid(krange,krange)
k=np.sqrt(k2d[0]**2+k2d[1]**2)
k=k.flatten()
fs=range(int(K/2))

tau=c*L
if gifDimensions == 1 and ntau>30: ntau=30

t0=0.
t_bound=tau*ntau
d=5
dt=dx/c/d

N=int((t_bound-t0)/dt)
print('N =',N)

#Set up variables
x=np.arange(0,L,dx)
y=np.arange(0,L,dy)
xx,yy = np.meshgrid(x,y)
u=np.zeros((K,K)) #f(t=n) #u[y,x]
U=np.zeros((K,K)) #f(t=n+1)
v=np.zeros((K,K))
V=np.zeros((K,K))
w=np.zeros((K,K))
W=np.zeros((K,K))
z=np.zeros((K,K))
Z=np.zeros((K,K))

ux=np.zeros((K,K))
uy=np.zeros((K,K))
t=np.linspace(t0,t_bound,N)

#Define funtions

#x-derivative
def Dx(u):
    for i in range(K):
        a=fft.fft(u[i,:])
        b=kx*a*1j
        ux[i,:]=np.real(fft.ifft(b))
    return ux

#y-derivative
def Dy(u):
    for i in range(K):
        ak=fft.fft(u[:,i])
        bk=ky*ak*1j
        uy[:,i]=np.real(fft.ifft(bk))
    return uy

#Differential equation
def f(a3,b3,c3,d3):
    a4=copy.deepcopy(b3)
    
    if equation=='WE': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3))
    if equation=='KG': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3)-m**2*a3)
    if equation=='SG': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3)-m**2*np.sin(np.pi/2*a3))
    if equation=='V4': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3)-m**2*a3-4*lam*a3**3)
    if equation=='uux': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3)-alpha*a3*(c3+d3))
    if equation=='uux2': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3)-beta*a3*(c3+d3)**2)
    if equation=='u2ux2': b4=copy.deepcopy(c**2*Dx(c3)+c**2*Dy(d3)-gamma*a3**2*(c3+d3)**2)
    
    c4=copy.deepcopy(Dx(b3))
    d4=copy.deepcopy(Dy(b3))
    return a4,b4,c4,d4

#Runge-Kutta 4
def RK4(a1,b1,c1,d1):    
    K1=f(a1,b1,c1,d1)
    K2=f(a1+dt/2*K1[0],b1+dt/2*K1[1],c1+dt/2*K1[2],d1+dt/2*K1[3])
    K3=f(a1+dt/2*K2[0],b1+dt/2*K2[1],c1+dt/2*K2[2],d1+dt/2*K2[3])
    K4=f(a1+dt*K3[0],b1+dt*K3[1],c1+dt*K3[2],d1+dt*K3[3])
    
    a2=copy.deepcopy(a1+dt/6*(K1[0]+2*K2[0]+2*K3[0]+K4[0]))
    b2=copy.deepcopy(b1+dt/6*(K1[1]+2*K2[1]+2*K3[1]+K4[1]))
    c2=copy.deepcopy(c1+dt/6*(K1[2]+2*K2[2]+2*K3[2]+K4[2]))
    d2=copy.deepcopy(d1+dt/6*(K1[3]+2*K2[3]+2*K3[3]+K4[3]))

    return a2,b2,c2,d2

#Calculation of the Power Spectrum
def PS(u):
    akk=np.abs(fft.fft2(u))
    a=akk.flatten()

    kbins = np.arange(0.5, K//2+1, 1.)

    Abins, _, _ = stats.binned_statistic(k, a,
                                         statistic = "mean",
                                         bins = kbins)
    Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

    Abins=np.flip(Abins)
    return(Abins)

#Set up initial conditions
u0 = pd.read_csv('u0.txt',header=None,sep=" ")
u0=u0.to_numpy()

u=u0
w=Dx(u)
z=Dy(u)

#Set up gif
if gifDimensions==2:
    fig = plt.figure()
    camera = Camera(fig)

if gifDimensions==3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
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
        if gifDimensions==2:
            plt.pcolormesh(x,y,u,cmap='inferno')
            camera.snap()     
            
        if gifDimensions==3:
            ax.plot_surface(xx, yy, u, cmap='inferno')
            camera.snap()
    
    #Calculate variables of the next time-step
    U,V,W,Z=RK4(copy.deepcopy(u),copy.deepcopy(v),copy.deepcopy(w),copy.deepcopy(z))
    
    #Break the loop if the simulation diverge
    if np.isnan(U[0,0])==True:
        break
    
    #Replace old variables with new variables
    u=copy.deepcopy(U)
    v=copy.deepcopy(V)
    w=copy.deepcopy(W)
    z=copy.deepcopy(Z)
    

print("Finished")

#Process gif
if gifDimensions==2 or gifDimensions==3:
    print("Processing gif...")
    if gifDimensions==3: print('Gif processing in 3D could take a couple of minutes, please wait')
    animation = camera.animate()
    animation.save('2D_Simulation.gif')
    plt.show()
    print("Gif finished")

#Save Power Spectrum
np.savetxt('PS_k.txt', fs, fmt='%d')
np.savetxt(r'PS_'+equation+'.txt', PS(u))

#Plot Power Spectrum
plt.plot(fs,PS(u),'--')
plt.xlabel(r"$k$")
plt.ylabel(r"$|a_k|$")
plt.yscale('log')
plt.tick_params(axis='both',direction='in',width=1,length=6)
plt.show()

