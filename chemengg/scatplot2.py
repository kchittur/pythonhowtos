from numpy import *
import pylab as p
from matplotlib import rc
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from scipy import stats, polyval, polyfit

R1a = 1
R2a = 10
KD1a = 0.1
KD2a = 5

def scatfA(x,a,b):
    return (a*x/(b + x))/x
def scatfB(x,a,b):
    return (a*x/(b + x))
def stline(x,slope,intcp):
    return x*slope + intcp

x = linspace(0.00001,200.0,4000)

lefty1a = scatfB(x,R1a,KD1a)
leftx1a = scatfA(x,R1a,KD1a)

lefty2a = scatfB(x,R2a,KD2a)
leftx2a = scatfA(x,R2a,KD2a)

ya = scatfB(x,R1a,KD1a) + scatfB(x,R2a,KD2a)
za = scatfA(x,R1a,KD1a) + scatfA(x,R2a,KD2a)

s1a, i1a, r1a, p1a, std_err = stats.linregress(ya[:3],za[:3])
s2a, i2a, r2a, p2a, std_err = stats.linregress(ya[-3:],za[-3:])

KD1Ca = -1.0/s1a
R1Ca = i1*KD1a
KD2Ca = -1.0/s2a
R2Ca = i2*KD2a

p.grid()
p.ylabel(r"\textbf{$\frac{LR}{L}$}",fontsize=16)
p.xlabel(r"\textbf{$[LR]$}")
p.xlim([0.0,8.0])
p.ylim([0.0,12.0])
p.plot(ya,za,lw=1)
p.plot(x,stline(x,s1a,i1a),lw=2,ls=':',color='r')
p.plot(x,stline(x,s2a,i2a),lw=2,ls=':',color='r')
p.title(r"Scatchard Plot for Two Receptors K$_{D1}$=%s~R$_{D1}$=%s~K$_{D2}$=%s~R$_{D2}$=%s" %(KD1,R1,KD2,R2),fontsize=12,color='b')  

p.savefig('concave-plot1.png')

p.clf()

