% noweb as the literate programming tool
% noweb enzyme-kinetics.w will extract the files you need
% and this make file can help
% name this enz.make

\documentclass{tufte-handout}
\usepackage{noweb}
\usepackage{fancyvrb}
\usepackage{minted}
\usepackage{listings}
\lstset{language=Python,prebreak=\raisebox{0ex}[0ex][0ex]
        {\ensuremath{\rhookswarrow}}}
\lstset{language=Python,postbreak=\raisebox{0ex}[0ex][0ex]
        {\ensuremath{\rcurvearrowse\space}}}
\lstset{language=Python,breaklines=true,breakatwhitespace=true}
\lstset{language=Python,numbers=left, numberstyle=\scriptsize,linewidth=6in,basicstyle=\scriptsize}
\usepackage{amsmath, amsthm, amssymb}
\usepackage{MnSymbol}
\usepackage{morefloats}

\usepackage{ifthen}

\newcommand{\executepython}{no}
\newcommand{\listpython}{no}

\usepackage{color}   %May be necessary if you want to color links
\usepackage{hyperref}
\hypersetup{
    colorlinks=true, %set true if you want colored links
    linktoc=all,     %set to all if you want both sections and subsections linked
    linkcolor=blue,  %choose some color if you want links to stand out
}

\usepackage{graphicx} % allow embedded images
\graphicspath{./}
%  \setkeys{Gin}{width=\linewidth,totalheight=\textheight,keepaspectratio}
%  \graphicspath{{graphics/}} % set of paths to search for images
\usepackage{amsmath}  % extended mathematics
\usepackage{booktabs} % book-quality tables
\usepackage{units}    % non-stacked fractions and better unit spacing
\usepackage{multicol} % multiple column layout facilities
\usepackage{lipsum}   % filler text
\usepackage{fancyvrb} % extended verbatim environments
  \fvset{fontsize=\normalsize}% default font size for fancy-verbatim environments

\usepackage{pythontex}
\setpythontexworkingdir{.}

\usepackage{enumitem,xcolor}

\usepackage[version=3]{mhchem}

\usepackage{multirow}
\usepackage{noweb}
\usepackage{eqnarray}

\usepackage{tikz}

\newcommand{\myred}[1]{\textcolor{red}{#1}}
\newcommand{\myblue}[1]{\textcolor{blue}{#1}}
\newcommand{\myeqn}[2]{\begin{equation*} #1 \tag{\bf {#2}} \label{#2} \end{equation*}}

\newcommand{\blue}[1]{{\textcolor{blue}{#1}}}
\newcommand{\red}[1]{{\textcolor{red}{#1}}}
\newcommand{\green}[1]{{\textcolor{green}{#1}}}

\hypersetup
{   pdfauthor = {Krishnan Chittur},
  pdftitle={Examples in Regression},
  colorlinks=TRUE,
  linkcolor=black,
  citecolor=blue,
  urlcolor=blue
}

\setlength{\parindent}{0pt}
\setlength{\parskip}{1.2ex}

\title{Analysis of Enzyme Kinetics}
\author{Krishnan Chittur \\ \url{http://webpages.uah.edu/~chitturk}}
\date{February 10, 2016}

\usepackage{comment}
\def\degF{$~^o$F~}
\def\degC{$~^o$C~}
\def\degR{$~^o$R~}

\begin{document}
\tableofcontents

\section{Abstract}
Code and explanation to create the McCabe Thiele Diagram for a binary system.

<<requiredlibraries>>=
'''
python script the McCabe Thiele binary distillation figure
author Krishnan K Chittur, chitturk@uah.edu
Use with caution ...
'''
import math
import sys
import matplotlib.pyplot as p
from scipy.optimize import fsolve 
import numpy as np
import os 

'''
fancy options - using TeX to annotate the figure(s) 
'''

from matplotlib import rc 
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import scipy.stats.mstats

@

\begin{equation}
P = \mbox{exp}\left[A - \frac{B}{T + C}\right]
\end{equation}

'''
define functions to calculate Pressure from Temperature and the
Antoine Coefficients
'''
<<PfromT>>=

def PfromT(A,B,C,T):
    return math.exp((A - B/(T + C)))

@

\begin{equation}
T = \frac{B}{A - \mbox{math}.\mbox{log}(P)}
\end{equation}

<<TfromP>>=

def TfromP(A,B,C,P):
    return B/(A - math.log(P)) - C

@

Calculate the xy plot - at a given Temperature

<<XYfromT>>=

'''
Generate an XY plot - you can calculate Py given Tx and generate the xy plot 
OR calculate Ty given Px and generate the xy plot
we begin with the Py from Tx
'''
def xyplotT(T,A1,B1,C1,A2,B2,C2):
    x1 = []
    y1 = []
    Ptotal = []
    npts = 100
    for n in range(npts):
        x = 0.0 + float(n)/float(npts)
        p1s = PfromT(A1,B1,C1,T)
        p2s = PfromT(A2,B2,C2,T)
        Pt = p1s*x + p2s*(1.0 - x)
        x1.append(x)
        y1.append(p1s*x/Pt)
        Ptotal.append(Pt)
    return x1, y1, Ptotal

@

Calculate the xy plot - at a given pressure

<<xyFromP>>=

'''
Now do the Ty given Px  
'''
def xyplotP(P,A1,B1,C1,A2,B2,C2):
    x1 = []
    y1 = []
    T = []
    npts = 100
    t1s = TfromP(A1,B1,C1,P)
    t2s = TfromP(A2,B2,C2,P)
    tmin = min(t1s,t2s)
    tmax = max(t1s,t2s)
    tint = (tmax - tmin)/float(npts)
    for n in range(npts):
        tadd = tmin + float(n)*tint
        p1s = PfromT(A1,B1,C1,tadd)
        p2s = PfromT(A2,B2,C2,tadd)
        x = (P - p2s)/(p1s - p2s)
        y = x*p1s/P
        x1.append(x)
        y1.append(y)
        T.append(tadd)
    return x1, y1, T

@

<<xyPwT>>=

'''
this function creates the xy plot, but also plots
a second y axis with the temperature 
'''
def xyplotPwT(P,A1,B1,C1,A2,B2,C2):
    import matplotlib.pyplot as plt
    fig1, ax1 = plt.subplots()
    x1 = []
    y1 = []
    T = []
    npts = 100
    t1s = TfromP(A1,B1,C1,P)
    t2s = TfromP(A2,B2,C2,P)
    tmin = min(t1s,t2s)
    tmax = max(t1s,t2s)
    tint = (tmax - tmin)/float(npts)
    for n in range(npts):
        tadd = tmin + float(n)*tint
        p1s = PfromT(A1,B1,C1,tadd)
        p2s = PfromT(A2,B2,C2,tadd)
        x = (P - p2s)/(p1s - p2s)
        y = x*p1s/P
        x1.append(x)
        y1.append(y)
        T.append(tadd)
    ax1.plot(x1,y1,'b-')
    ax1.set_xlabel(' Mole Fraction (x) ')
    ax1.set_ylabel(' Mole Fraction (y) ')
    ax2 = ax1.twinx()
    ax2.set_ylabel(' Temperature ')
    ax2.plot(x1,T,color='r')
    plt.show() 

@

At different times, we have to determine the intersection points of
two straight lines - this function does that

<<intersection>>=

'''
calculate the intersection of two straight lines 
(x1,y1) to (x2,y2) and (x3,y3) to (x4,y4) 
we need it when drawing the mccabe thiele diagram
'''
def intersection(x1,y1,x2,y2,x3,y3,x4,y4):  
 one = (x1*y2 - y1*x2)*(x3 - x4) 
 two = (x1 - x2)*(x3*y4 - y3*x4) 
 three = (x1*y2 - y1*x2)*(y3 - y4) 
 four = (y1 - y2)*(x3*y4 - y3*x4) 
 five = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)
 Px = (one - two)/five
 Py = (three - four)/five  
 return Px, Py  

@

Now create a function to do the McCabe Thiele Calculations

<<mccabethiele>>=

'''
assemble everything to do the mccabe thiele
'''
def mccabe(P,xf,xd,xb,q,Rf):
    '''
    start with creating an XY plot at a given Pressure (so the Temperature will vary)
    '''
    x = xyplotP(P,A1,B1,C1,A2,B2,C2)[0]
    y = xyplotP(P,A1,B1,C1,A2,B2,C2)[1]
    T = xyplotP(P,A1,B1,C1,A2,B2,C2)[2]
    '''
    now calculate an average alpha
    '''
    myalpha = []
    for i in range(len(x)):
        myalpha.append(y[i]*(1.0 - x[i])/((1.0 - y[i])*x[i]))
    myalphaavg = np.average(myalpha)   
    '''
    here is some code to calculate the geometric average ...
#    geomavg = scipy.stats.mstats.gmean(myalpha)
#    print myalphaavg, geomavg
    '''    
    p.grid()

# make sure the x and y axes are of the same size

    p.axes().set_aspect('equal')
#   plot the 45 degree line    
    x1 = np.linspace(0.0,1.0,100)
    y1 = x1
    p.xlim([0.0,1.0])
    p.ylim([0.0,1.0])
    p.xlabel('Mole Fraction (x)')
    p.ylabel('Mole Fraction (y)')
    p.title('XY at P=%s' %(P))
    '''
    plot the xy diagram - save the figures along the way so the animation 
    can be created at the end
    '''
    p.plot(x, y, 'b', label = 'x y', lw = 1)
    p.title(r"McCabe Thiele for a Binary System ") 
    p.text(0.1,0.95,r"$x_D$ = %s, $x_B$ = %s, $x_F$ = %s, q = %s, $R_f$ = %s " %(xd, xb, xf, q, Rf),fontsize=8) 
    p.savefig('mct1.png')
    p.title(r"XY Plot using Antoine''s Equation and Raoult''s Law for %s " %(systemname),fontsize=8,color='b')    
    '''
    plot the VLE using average alpha
    '''
    y2 = []
    for i in range(len(x)):
        y2.append(myalphaavg*x[i]/(1.0 + (myalphaavg - 1.0)*x[i]))

    p.plot(x,y2,'r')
    p.savefig('mct2.png')
    
    p.title(r"The VLE XY Plot overlayed by $y=\frac{\alpha x}{1.0 + (\alpha - 1.0)x}$",fontsize=12,color='b')
    p.savefig('mct3.png')
    
    p.plot(x1, y1, 'b', label = '', linestyle = 'dashed', lw = 2)
    p.title(r"Plotting the 45 degree line",fontsize=12)
    p.savefig('mct4.png')
    
    '''
    locate xf, xb and xd
    '''
    markers_on = [xd*100.0]
#    markers_on = xd*100.0
    p.plot(x1, y1, 'g-')
    p.title(r' Locating the Distillate Composition $x_D$ = %s'%(xd),fontsize=12)
    p.plot(x1, y1, 'ro', markevery=markers_on)
#    p.plot(x1[markers_on], y1[markers_on], 'ro')
    p.savefig('mct5.png')
    
    markers_on = [xf*100.0]
#    markers_on = xf*100.0
    p.title(r'Locating the Feed Composition $x_F$ = %s'%(xf),fontsize=12)
    p.plot(x1,y1,'bo', markevery=markers_on)
#    p.plot(x1[markers_on], y1[markers_on], 'bo')
    p.savefig('mct6.png')
    
    markers_on = [xb*100.0]
    p.title(r'Locating the bottoms composition $x_B$ = %s'%(xb),fontsize=12)
    p.plot(x1,y1,'ro', markevery=markers_on)
#    p.plot(x1[markers_on], y1[markers_on],'ro')
    p.savefig('mct7.png')
    '''
    now plot the q line, calculate the slope using the q value
    '''
    qslope = q/(q - 1.0)
    '''
    now we determine where the q line intersects the xy plot 
    '''
    def fun1(x):
     return myalphaavg*x/(1.0 + (myalphaavg - 1.0)*x)
    def fun2(x):
     return qslope*x + -xf/(q-1)
    def findIntersection(fun1,fun2,x0):
     return fsolve(lambda x : fun1(x) - fun2(x),x0)

    qintx = findIntersection(fun1,fun2,0.0)
    qinty = fun1(qintx) 
    p.title(r'Now plotting the q line for q = %s, q slope = $\frac{q}{q-1}$ = %s'%(q,qslope),fontsize=10)
    p.plot([xf,qintx],[xf,qinty])
    p.savefig('mct8.png')
    p.title(r'Joining the intersection of the q line with the XY plot to $x_D$= %s'%(xd),fontsize=10)
    p.plot([qintx,xd],[qinty,xd])
    p.savefig('mct9.png')
    rminslope = (xd - qinty)/(xd - qintx)
    rminint = qinty - rminslope*qintx
    Rmin = xd/rminint - 1.0 
    Ractual = Rmin*Rf
    actint = xd/(Ractual + 1) 
    p.title(r'Extending the $R_{min}$ line to Y axis to get intercept $\frac{x_D}{R_{min}+1}$ = %s'%(rminint),fontsize=10)
    p.plot([0.0,qintx],[rminint,qinty])
    p.savefig('mct10.png')
    p.title(r'Multiply $R_{min}$=%s by $Rf=$%s to get new Intercept $\frac{x_D}{R+1}$ = %s'%(Rmin,Rf,actint),fontsize=10)
    p.plot([0.0,xd],[actint,xd])
    p.savefig('mct11.png')
    p.plot([0.0],[rminint],'rH')
    p.savefig('mct12.png')
    p.plot([0.0],[actint],'rD')
    p.savefig('mct13.png')
    p.plot([0.0,xd],[actint,xd])
    p.savefig('mct14.png')
    '''
    determine intersection 
    '''
    bx = intersection(0.0,actint,xd,xd,xf,xf,qintx,qinty)[0]
    by = intersection(0.0,actint,xd,xd,xf,xf,qintx,qinty)[1]
    rbx = round(bx,2)
    rby = round(by,2)
    p.title(r'Join $x_B$=%s to intersection point $x$=%s, $y=$%s'%(xb,rbx,rby),fontsize=10)
    p.plot([xb,bx],[xb,by])
    p.savefig('mct15.png')
    p.plot([bx],[by],'ro')
    p.savefig('mct16.png')
#    p.grid()
    roplineslope = (xd - actint)/xd
    soplineslope = (by - xb)/(bx - xb)
    ropline = roplineslope*x + actint
    sopline = soplineslope*x + xb*(1.0 - soplineslope)
    nr = 0 
    ns = 0
    '''
    start the stepping for the stages
    calculate the x intersection point of the first horizontal line at xd
    start the stepping for the stage, at xd 
'''
    myn = 17
    p.title(r'Start calculating the number of Stages')
    nsnames = 'ns' 
    ssnames = 'ss'
    myx = xd
    myy = xd
    while ((myx >= xb) and (myy >= xb)):
        myx1 = myy/(myalphaavg - (myalphaavg - 1.0)*myy)
        myy1 = myy
        p.plot([myx,myx1],[myy,myy1],ls=':',lw=2,color='b')
        p.savefig('mct'+str(myn)+'.png')
#        myn = myn + 1
        myx = myy/(myalphaavg - (myalphaavg - 1.0)*myy)
        if ((myx >= bx) and (myy > by)):
            p.title(r" Rectification Section ")
            myx = myx1
            myy = myy1 
            myx1 = myx
            myy1 = myalphaavg*myx/(1.0 + (myalphaavg - 1.0)*myx)
            p.plot([myx,myx1],[myy,myy1],ls=':',lw=2,color='b')
            myn = myn + 1
            p.savefig('mct'+str(myn)+'.png')
            
#            p.axhline(y=myy,ls=':',lw=1,color='r')            
            myy = roplineslope*myx + actint
            myy1 = myalphaavg*myx/(1.0 + (myalphaavg - 1.0)*myx)
            p.plot([myx,myx],[myy1,myy],ls=':',lw=2,color='b')
            myn = myn + 1
            p.savefig('mct'+str(myn)+'.png')
#            p.axvline(x=myx,ls=':',lw=1,color='r')
            nr = nr + 1
            p.title(r"\textit{Rectification stages} = %s" %(nr),fontsize=12,color='b')
            nsn = str(nsnames)+str(nr)+'.png'
            p.savefig(nsn)
        else:
        
#            p.axhline(y=myy,ls=':',lw=1,color='r') 
#            p.title(r" Start Calculating the Stripping Stages ")
            p.title(r"\textit{Stripping section}")
            myx = myy/(myalphaavg - (myalphaavg - 1.0)*myy)
            x1 = myx
            y1 = myy
            myy = soplineslope*myx + xb*(1.0 - soplineslope)
            y2 = myy
#            p.axvline(x=myx,ls=':',lw=1,color='b')
            p.plot([x1,x1],[y2,y1],ls=':', lw=2,color='b')
            myn = myn + 1
            p.savefig('mct'+str(myn)+'.png')
            ns = ns + 1 
            p.title(r"$N_{rect}$ = %s $N_{stripping}$ = %s" %(nr, ns))
            ssn = str(ssnames)+str(ns)+'.png'
            p.savefig(ssn)

    minstages = 0 
    '''
    start the stepping for the stages
    calculate the x intersection point of the first horizontal line at xd
    '''
    p.title(r" Calculate the minimum number of stages ")
    myx = xd
    myy = xd
    while ((myx >= xb) and (myy >= xb)):
        myx = myy/(myalphaavg - (myalphaavg - 1.0)*myy)
        myn = myn + 1 
        
        p.plot([myx,myy],[myy,myy],ls='-.',lw=3,color='r')    
        p.savefig('mct'+str(myn)+'.png')
#        p.axhline(y=myy,ls=':',lw=2,color='r')
        myn = myn + 1
        myy = myalphaavg*myx/(1.0 + (myalphaavg - 1.0)*myx)
        p.plot([myx,myx],[myx,myy],ls='-.',lw=3,color='r')
        p.savefig('mct'+str(myn)+'.png')
#        p.axvline(x=myx,ls=':',lw=1,color='r')
        minstages = minstages + 1
        myy = myx



#    print "myn now is ", myn 
    p.title('Rectification Stages = %s, Stripping Stages = %s' %(nr,ns))
    p.text(0.1,0.9,r"$x_d$ = %s, $x_b$ = %s, $x_f$ = %s" %(xd,xb,xf))
    p.text(0.1,0.8,r"q = %s, $R_f$ = %s" %(q,Rf))
    p.text(0.3,0.1,r"$N_{min} = %s N_{rect}$ = %s, $N_{strip}$ = %s " %(minstages,nr,ns))
    p.text(0.4,0.3,r"$\alpha_{avg}$ = %s " %(round(myalphaavg,4)))
    p.text(0.3,0.05,r"$R_{min}$ = %s, $R_{actual}$ = %s " %(Rmin, Ractual))
    Tlow = round(min(T),1)
    Thigh = round(max(T),1)
    p.text(0.3,0.2,r"XY at P = %s, $T_{low}$ = %s, $T_{high}$ = %s " %(P,Tlow,Thigh))
    p.text(0.45,0.4,r"System Name = %s " %(systemname),fontsize=10)
    p.text(-0.2,1.3*rminint,r"$\frac{%s}{%s + 1}$" %(xd,round(rminint,3)))
    myn = myn + 1
#    print "myn is now ", myn
    p.savefig('mct'+str(myn)+'.png')
    p.annotate("", xy=(0.0, 1.02*rminint),  xycoords='data',
                xytext=(-50, 50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc,angleA=10,armA=5,rad=5"),
                )
    
#    p.savefig('mct'+str(myn)+'.png')
#    print "nr = ", nr, "ns = ", ns
    p.show()
    p.clf()
    return

@

<<binarysystem>>=

#http://stackoverflow.com/questions/8409095/matplotlib-set-markers-for-individual-points-on-a-line

'''
now create some plots 
these values are for benzene (1), and toluene (2) from Van Ness/7th edition
P in kPa, T in degrees C
'''
A1 = 13.7819
B1 = 2726.81
C1 = 217.572
A2 = 13.9320
B2 = 3056.96
C2 = 217.625
P = 101.300    
systemname = 'Benzene-Toluene'

@

<<specifications>>=

xf = 0.40
xd = 0.95
xb = 0.05
q = 0.99
Rf = 1.01

@

<<specifications1>>=

xf = 0.40
xd = 0.95
xb = 0.05
q = 0.99
Rf = 1.21

@



<<compileimages>>=

from PIL import Image

import glob
import os
mymctpng = []
for files in glob.glob("mct*.png"):
    mymctpng.append(files)
mynspng = []
for files in glob.glob("ns*.png"):
    mynspng.append(files)
mysspng = []
for files in glob.glob("ss*.png"):
    mysspng.append(files)

print (mymctpng, len(mymctpng))

from natsort import natsorted

mymctpngsorted = natsorted(mymctpng,reverse=False)
mynspngsorted = natsorted(mynspng,reverse=False)
mysspngsorted = natsorted(mysspng,reverse=False)

import re


from os.path import basename

mygif = []

for i in range(len(mymctpngsorted)):
    im = Image.open(mymctpngsorted[i])
    alpha = im.split()[3]
    im = im.convert('RGB').convert('P', palette=Image.ADAPTIVE, colors=255)
    mask = Image.eval(alpha, lambda a: 255 if a <=128 else 0)
    im.paste(255, mask)
    outname = os.path.splitext(mymctpngsorted[i])[0]+'.gif'
    print ("Converted ", mymctpngsorted[i], " to ", outname)
    mygif.append(outname)
    im.save(outname, transparency=255)

for i in range(len(mynspngsorted)):
    im = Image.open(mynspngsorted[i])
    alpha = im.split()[3]
    im = im.convert('RGB').convert('P', palette=Image.ADAPTIVE, colors=255)
    mask = Image.eval(alpha, lambda a: 255 if a <=128 else 0)
    im.paste(255, mask)
    outname = os.path.splitext(mynspngsorted[i])[0]+'.gif'
    print ("Converted ", mynspngsorted[i], " to ", outname)
    mygif.append(outname)
    im.save(outname, transparency=255)

for i in range(len(mysspngsorted)):
    im = Image.open(mysspngsorted[i])
    alpha = im.split()[3]
    im = im.convert('RGB').convert('P', palette=Image.ADAPTIVE, colors=255)
    mask = Image.eval(alpha, lambda a: 255 if a <=128 else 0)
    im.paste(255, mask)
    outname = os.path.splitext(mysspngsorted[i])[0]+'.gif'
    print ("Converted ", mysspngsorted[i], " to ", outname)
    mygif.append(outname)
    im.save(outname, transparency=255)

def animate(fname, max_gen):
   string = ""
   mgen = max_gen 
   for i in range(1, mgen):
      if os.stat(str(fname) + str(i) + ".gif"):
         string += str(fname) + str(i) + ".gif "
   os.system("convert" + ' ' + '-delay 200 -loop 0 ' + string + "animation"+str(fname)+".gif")
   print ("Finished animate")

#print len(mymctpngsorted)

animate('mct',len(mymctpngsorted)) 
#animate('ns',len(mynspngsorted))
#animate('ss',len(mysspngsorted))

os.system("mv animationmct.gif mccabe-thiele-example5.gif") 
os.system("ffmpeg -f image2 -r 1/5 -i mct%01d.png -vcodec mpeg4 -y mccabe-thiele-example5.mp4")
 
@

<<mccabethieleplot.py>>=
<<requiredlibraries>>
<<PfromT>>
<<TfromP>>
<<XYfromT>>
<<xyFromP>>
<<xyPwT>>
<<intersection>>
<<mccabethiele>>
<<binarysystem>>
<<specifications1>>
mccabe(P,xf,xd,xb,q,Rf)
<<compileimages>>

#xyplotPwT(P,A1,B1,C1,A2,B2,C2)

'''



'''
