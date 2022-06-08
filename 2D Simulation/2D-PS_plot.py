import matplotlib.pyplot as plt
import pandas as pd


plt.rc( 'font', size=18, family="DejaVu Sans" )
plt.rc( 'text', usetex=True) #per format LaTex

# Run only if all the txt files are created with the same K

k = pd.read_csv('PS_k.txt',header=None)
WE = pd.read_csv('PS_WE.txt',header=None)
KG = pd.read_csv('PS_KG.txt',header=None)
sineG = pd.read_csv('PS_SG.txt',header=None)
V4 = pd.read_csv('PS_V4.txt',header=None)
uux = pd.read_csv('PS_uux.txt',header=None)
uux2 = pd.read_csv('PS_uux2.txt',header=None)
u2ux2 = pd.read_csv('PS_u2ux2.txt',header=None)


plt.plot(WE,'k:',label="Wave Eq")
plt.plot(KG,'--',label='KG Eq')
plt.plot(sineG,'--',label='sineG Eq')
plt.plot(V4,'--',label='$V=\lambda u^4$')
plt.plot(uux,'--',label=r'$u\nabla u$')
plt.plot(uux2,'--',label=r'$u(\nabla u)^2$')
plt.plot(u2ux2,'--',label=r'$u^2(\nabla u)^2$')
plt.xlabel(r"$k$")
plt.ylabel(r"$|a_k|$")
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.yscale('log')
plt.tick_params(axis='both',direction='in',width=1,length=6)
plt.show()

