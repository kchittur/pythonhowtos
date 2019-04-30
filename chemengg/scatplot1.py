from numpy import *
import pylab as p
from matplotlib import rc
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from scipy import stats, polyval, polyfit

R1 = 1
R2 = 10
KD1 = 0.1
KD2 = 10.0

def scatfA(x,a,b):
    return (a*x/(b + x))/x
def scatfB(x,a,b):
    return (a*x/(b + x))
def stline(x,slope,intcp):
    return x*slope + intcp

x = linspace(0.00001,200.0,4000)

lefty1 = scatfB(x,R1,KD1)
leftx1 = scatfA(x,R1,KD1)

lefty2 = scatfB(x,R2,KD2)
leftx2 = scatfA(x,R2,KD2)

y = scatfB(x,R1,KD1) + scatfB(x,R2,KD2)
z = scatfA(x,R1,KD1) + scatfA(x,R2,KD2)

s1, i1, r1, p1, std_err = stats.linregress(y[:3],z[:3])
s2, i2, r2, p2, std_err = stats.linregress(y[-3:],z[-3:])

KD1C = -1.0/s1
R1C = i1*KD1
KD2C = -1.0/s2
R2C = i2*KD2

p.grid()
p.ylabel(r"\textbf{$\frac{LR}{L}$}",fontsize=16)
p.xlabel(r"\textbf{$[LR]$}")
p.xlim([0.0,8.0])
p.ylim([0.0,12.0])
#p.plot(leftx1,lefty1,lw=1,color='b')
#p.plot(leftx2,lefty2,lw=1,ls=':',color='b')
p.plot(y,z,lw=1)
p.plot(x,stline(x,s1,i1),lw=2,ls=':',color='r')
p.plot(x,stline(x,s2,i2),lw=2,ls=':',color='r')
p.title(r"Scatchard Plot for Two Receptors K$_{D1}$=%s~R$_{D1}$=%s~K$_{D2}$=%s~R$_{D2}$=%s" %(KD1,R1,KD2,R2),fontsize=12,color='b')  

p.show()

#p.savefig('concave-plot.png')

p.clf()

#print R1, R2, KD1, KD2
#print R1C, R2C, KD1C, KD2C
