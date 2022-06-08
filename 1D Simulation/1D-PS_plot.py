import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import fftpack as fft

# Run only if all the txt files are created with the same K

k = fft.fftshift(pd.read_csv('PS_k.txt',header=None)[0])
WE = fft.fftshift(pd.read_csv('PS_WE.txt',header=None)[0])
KG = fft.fftshift(pd.read_csv('PS_KG.txt',header=None)[0])
sineG = fft.fftshift(pd.read_csv('PS_SG.txt',header=None)[0])
V4 = fft.fftshift(pd.read_csv('PS_V4.txt',header=None)[0])
uux = fft.fftshift(pd.read_csv('PS_uux.txt',header=None)[0])
uux2 = fft.fftshift(pd.read_csv('PS_uux2.txt',header=None)[0])
u2ux2 = fft.fftshift(pd.read_csv('PS_u2ux2.txt',header=None)[0])

k=k[-int(len(k)/2):]
WE=WE[-int(len(WE)/2):]
KG=KG[-int(len(KG)/2):]
sineG=sineG[-int(len(sineG)/2):]
V4=V4[-int(len(V4)/2):]
uux=uux[-int(len(uux)/2):]
uux2=uux2[-int(len(uux2)/2):]
u2ux2=u2ux2[-int(len(u2ux2)/2):]


plt.plot(k,WE,'k:',label="Wave Eq")
plt.plot(k,KG,'--',label='KG Eq')
plt.plot(k,sineG,'--',label='sineG Eq')
plt.plot(k,V4,'--',label='V=$\lambda·u^4$')
plt.plot(k,uux,'--',label='$u·\partial_x u$')
plt.plot(k,uux2,'--',label='$u·(\partial_x u)^2$')
plt.plot(k,u2ux2,'--',label='$u^2·(\partial_x u)^2$')
plt.xlabel(r"$k$")
plt.ylabel(r"$|a_k|$")
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.yscale('log')
plt.tick_params(axis='both',direction='in',width=1,length=6)
plt.show()



