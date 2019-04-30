% https://www.polypompholyx.com/2013/05/misconceptions-mmm/ (mm kinetics not valid inside cell)
% https://www.ncbi.nlm.nih.gov/books/NBK22430/

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
The enzymatic conversion of para nitrophenol phosphate (pNPP)
to para nitrophenol (pNP) with an alkaline phosphatase has been well studied.
This reaction is also frequently used to explain and use the Michaelis
Menten rate equation for analysis of the conversion of the pNPP (colorless)
to pNP (yellow color).  What we have attempted to do here is to start with 
a detailed model of the enzymatic reaction and approximations that lead
to the Michaelis Menten equation.  We provide computer programs (in python)
that will allow students to examine assumptions involved in deriving the
MM equation, fit the models to experimental data and develop a better
understanding of enzyme kinetics.
\subsection{The document and the contents}
I have used what is termed as {\em literal programming} - where I have 
mixed documentation, code, sample data files together and used a program called
{\em noweb} to create/extract the documentation, the code and the data files.
There is also a {\em make} file that can simplify the process of compiling the
documentation.

\subsection{The Make file}
A {\em make} file contains instructions on how to assemble the code, documentation
(or whatever you are trying to do).  In this example, we have {\em code} (in python),
documentation in \LaTeX - we have also used 
pythontex which allows for the embedding of python within \LaTeX and display
results from the execution of the python code.

\subsection{How to extract and compile}
You will need \href{https://www.cs.tufts.edu/~nr/noweb/}{noweb} if you begin with this
file - the command {\em noweb analyze-enzyme-kinetics.w} will extract the
\TeX file, the python codes and the make file.  You can then run
{\em make -f analyze.make} to {\bf make} the documentation and run the code.
Sample data files are also included here.

<<enz.make>>=
all: enzyme-kinetics.pdf 

enzyme-kinetics.tex: enzyme-kinetics.w
	noweb enzyme-kinetics.w

directfit.png: directfitPvsTime.py
	python2 directfitPvsTime.py

directfitS0fixed.png: directfitPvsTimeSfixed.py
	python2 directfitPvsTimeSfixed.py

lvbplot.png: initialrates.py
	python2 initialrates.py

enzyme-kinetics.pdf: enzyme-kinetics.tex directfit.png directfitS0fixed.png lvbplot.png
	pdflatex -shell-escape enzyme-kinetics.tex
	pythontex --interpreter python:python2 enzyme-kinetics.tex
	pdflatex -shell-escape enzyme-kinetics.tex
	pdflatex -shell-escape enzyme-kinetics.tex

clean:
	rm -r -f *.aux *.toc enzyme-kinetics.pdf enzyme-kinetics.tex _minted*
	rm -r -f pythontex-files-enzyme-kinetics
	rm -f enzyme-kinetics.pytxcode

@

\subsection{Required Libraries}
A collection of python libraries that many of the scripts will need 

<<requiredlibraries>>=
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import scipy
from scipy import *
from scipy.integrate import odeint
from scipy import integrate
from scipy.optimize import fmin
from scipy.optimize import curve_fit
from scipy import stats, polyval, polyfit  
# matplotlib allows the use of TeX for equations in plots
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import os.path

@


\subsection{Reading blank lines}
often we do not know how many lines a particular data file will have - here we 
read a file line by line till it reaches a blank or the end

<<readnonblanklines>>=
def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

@


\subsection{Convert absorbance to concentration}


<<absorbance2concentration>>=
def abs2conc(infile):
    inf = open(infile,'r')
    t = []
    a = []
    nl = 0
    for line in nonblank_lines(inf):
        nl = nl + 1 
        thisline = line.split();
        t.append(float(thisline[0]))
        a.append(float(thisline[1])/18800.0)
    return t, a

@



\section{Some data files}
\subsection{Data files used in these calculations} Sample files are here

<<zs1e1.dat>>=
0       0.008
20      0.124
40      0.255
60      0.372
80      0.475
100     0.57
120     0.655
140     0.742
160     0.835

@

<<zs1e2.dat>>=
0       0.007
20      0.065
40      0.124
60      0.182
80      0.242
100     0.304
120     0.365
140     0.424
160     0.476
180     0.526
200     0.58
220     0.636
240     0.691
260     0.745
280     0.796
300     0.845

@

<<zs1e3.dat>>=
0       0.007
20      0.038
60      0.099
80      0.13
100     0.159
120     0.192
140     0.222
160     0.25
180     0.278
200     0.308
220     0.338
240     0.365
260     0.392
280     0.42
300     0.45

@

<<zs1e4.dat>>=
0       0.007
20      0.02
40      0.033
60      0.048
80      0.065
100     0.082
120     0.099
140     0.118
160     0.136
180     0.151
200     0.167
220     0.179
240     0.189
260     0.203
280     0.219
300     0.238

@

<<zs2e1.dat>>=
0       0.004
20      0.123
40      0.241
60      0.35
80      0.451
100     0.552
120     0.64
140     0.734
160     0.825

@

<<zs2e2.dat>>=
0       0.004
20      0.069
40      0.131
60      0.19
80      0.247
100     0.301
120     0.385
140     0.418
160     0.476
180     0.53
200     0.582
220     0.631
240     0.678
260     0.723
280     0.792
300     0.815

@

<<zs2e3.dat>>=
0       0.003
20      0.032
40      0.045
60      0.09
80      0.122
100     0.153
120     0.185
140     0.212
160     0.24
180     0.27
200     0.302
220     0.334
240     0.366
260     0.396
280     0.426
300     0.454

@

<<zs2e4.dat>>=
0       0.003
20      0.021
40      0.037
60      0.052
80      0.067
100     0.081
120     0.098
140     0.115
160     0.131
180     0.147
200     0.163
220     0.179
240     0.195
260     0.211
280     0.227
300     0.244

@

<<zs3e1.dat>>=
0       0.001
20      0.12
40      0.232
60      0.327
80      0.414
100     0.492
120     0.563
140     0.628
160     0.69
180     0.742
200     0.787
220     0.83
240     0.865

@

<<zs3e2.dat>>=
0       0.002
20      0.066
40      0.126
60      0.18
80      0.23
100     0.278
120     0.327
140     0.38
160     0.424
180     0.466
200     0.506
220     0.54
240     0.574
260     0.608
300     0.674

@

<<zs3e3.dat>>=
0       0.002
20      0.032
40      0.06
60      0.089
80      0.119
100     0.149
120     0.177
140     0.204
160     0.23
180     0.255
200     0.28
220     0.307
240     0.333
260     0.357
280     0.38
300     0.402

@

<<zs3e4.dat>>=
0       0.001
20      0.016
40      0.032
60      0.047
80      0.062
100     0.076
120     0.09
140     0.104
160     0.118
180     0.132
200     0.147
220     0.162
240     0.176
260     0.19
280     0.203
300     0.216

@

<<zs4e1.dat>>=
0       0.003
20      0.11
40      0.194
60      0.262
80      0.314
100     0.358
120     0.392
140     0.418
160     0.44
180     0.458
200     0.472
220     0.486
240     0.494
260     0.502
280     0.508
300     0.51

@

<<zs4e2.dat>>=
0       0.002
20      0.05
40      0.096
60      0.138
80      0.176
100     0.211
120     0.243
140     0.27
160     0.295
180     0.317
200     0.339
220     0.358
240     0.375
260     0.39
280     0.404
300     0.416

@

<<zs4e3.dat>>=
0       0.003
20      0.03
40      0.059
60      0.086
80      0.112
100     0.137
120     0.16
140     0.182
160     0.202
180     0.22
200     0.239
220     0.255
240     0.271
260     0.285
280     0.299
300     0.313

@

<<zs4e4.dat>>=
0       0.001
20      0.016
40      0.029
60      0.042
80      0.056
100     0.07
120     0.085
140     0.098
160     0.11
180     0.123
200     0.135
220     0.146
240     0.156
260     0.165
280     0.175
300     0.186

@

<<zs5e1.dat>>=
0       0.001
20      0.085
40      0.144
60      0.183
80      0.21
100     0.228
120     0.241
140     0.25
160     0.254
180     0.261
200     0.264
220     0.266
240     0.267
260     0.268
280     0.268
300     0.269

@

<<zs5e2.dat>>=
0       0.001
40      0.09
60      0.123
80      0.148
100     0.169
120     0.188
140     0.203
160     0.216
180     0.226
200     0.234
220     0.24
240     0.246
260     0.251
280     0.254
300     0.257

@

<<zs5e3.dat>>=
0       0.001
20      0.029
40      0.054
60      0.074
80      0.093
100     0.11
120     0.126
140     0.143
180     0.17
200     0.18
220     0.189
240     0.197
260     0.204
280     0.211
300     0.217

@

<<zs5e4.dat>>=
0       0
20      0.014
40      0.029
60      0.042
80      0.053
100     0.064
120     0.075
140     0.088
160     0.099
180     0.109
200     0.118
220     0.125
240     0.132
260     0.14
280     0.147
300     0.154

@

<<zs6e1.dat>>=
0       0.001
20      0.071
40      0.104
60      0.12
80      0.128
100     0.132
120     0.134
140     0.136
160     0.136
180     0.136

@

<<zs6e2.dat>>=
0       0.002
20      0.045
40      0.074
60      0.094
80      0.107
100     0.116
120     0.122
140     0.126
160     0.129
180     0.13
200     0.132
220     0.132
240     0.133
260     0.133
280     0.133
300     0.134

<<zs6e3.dat>>=
0       0.001
20      0.023
40      0.041
60      0.055
80      0.068
100     0.078
120     0.087
140     0.096
160     0.103
180     0.108
200     0.112
220     0.115
240     0.118
260     0.12
280     0.123
300     0.125

<<zs6e4.dat>>=
0       0.001
20      0.011
40      0.023
60      0.033
80      0.041
100     0.049
120     0.056
140     0.061
160     0.064
180     0.069
200     0.075
220     0.081
240     0.083
260     0.087
280     0.089
300     0.092

@

<<zconcentration1.dat>>=
zs1e1.dat 0.0003 15.36e-9
zs2e1.dat 0.00015 15.36e-9
zs3e1.dat 0.000075 15.36e-9
zs4e1.dat 0.00003 15.36e-9
zs5e1.dat 0.000015 15.36e-9
zs6e1.dat 0.0000075 15.36e-9

@

<<zconcentration2.dat>>=
zs1e2.dat 0.0003 7.68e-9
zs2e2.dat 0.00015 7.68e-9
zs3e2.dat 0.000075 7.68e-9
zs4e2.dat 0.00003 7.68e-9
zs5e2.dat 0.000015 7.68e-9
zs6e2.dat 0.0000075 7.68e-9

@

<<zconcentration3.dat>>=
zs1e3.dat 0.0003 3.84e-9
zs2e3.dat 0.00015 3.84e-9
zs3e3.dat 0.000075 3.84e-9
zs4e3.dat 0.00003 3.84e-9
zs5e3.dat 0.000015 3.84e-9
zs6e3.dat 0.0000075 3.84e-9

@

<<zconcentration4.dat>>=
zs1e4.dat 0.0003 1.92e-9
zs2e4.dat 0.00015 1.92e-9
zs3e4.dat 0.000075 1.92e-9
zs4e4.dat 0.00003 1.92e-9
zs5e4.dat 0.000015 1.92e-9
zs6e4.dat 0.0000075 1.92e-9

@

\section{This document}
Our objective is to describe as completely as possible an experiment and calculations
related to the experiment involving an enzyme substrate reaction.  The enzyme is 
alkaline phosphatase and the substrate is para nitrophenol phosphate.
This system is used in high schools and colleges to explain several core ideas
in enzyme reactions.  We will start with the conventional approaches, explanations 
followed by extensive and detailed analysis of the equations and analysis of data.

\section{The Enzyme Mechanism, as normally presented}
The mechanism proposed for the conversion of a substrate $S$
to a product, $P$ is shown in ~\ref{conventional}.  

\begin{equation}
\ce{E + S <=>T[\cf{k_1}][\cf{k_{-1}}] E--S ->T[\cf{k_{2}}] E + P }
\label{conventional}
\end{equation}

The substrate is postulated to bind
to the enzyme (E) forming an enzyme 
substrate complex (E--S) which then breaks up into the original enzyme
and the product (P)

A widely accepted approach to arrive at a {\em rate equation} from the
mechanism shown as ~\ref{conventional}
uses either what is called a quasi steady state approximation (QSSA) or 
an equilibrium assumption (EA).  In the QSSA, the assumption is made that 
the enzyme substrate complex concentration does not change appreciably
with time ($\frac{d[ES]}{dt} ~\approx 0.0$)  In the
EA, the assumption is that the binding of the enzyme to the free substrate 
happens rapidly and an equilibrium is established 
with the enzyme substrate complex ($\frac{[E][S]}{[ES]} = K_{eq}$)
Details
of these approximations and derivations 
can be found in many places.  The final form of the rate 
equation is called the Michaelis Menten enzymatic rate equation 
and is shown ~\ref{michaelismentenrate}

\begin{equation}
\mbox{rate} = V = \frac{dP}{dt} = -\frac{dS}{dt} = \frac{V_m[S]}{[S] + K_m}
\label{michaelismentenrate}
\end{equation}

where $V_m = k_2 E_0$, $K_m = \frac{k_{-1} + k_2}{k_1}$ (for QSSA) 
and $K_m = \frac{k_{-1}}{k_1}$ (for EA).  
$E_0$ is the initial concentration of the enzyme and $S$ the concentration
of the substrate (which is a function of time).  

Typical enzyme kinetics experiments are designed to use 
the equation~\ref{michaelismentenrate} to calculate the
$V_m$ and $K_m$.  This is done by 
using what is called the
Lineweaver Burke Plot (also called the double reciprocal plot)

\begin{equation}
\frac{1}{\mbox{rate}} = \frac{1}{V} = \frac{1}{\frac{dP}{dt}} = \frac{K_m}{V_m}\frac{1}{S} + \frac{1}{V_m}
\label{lwbequation}
\end{equation}

We can see from equation~\ref{lwbequation} that a plot of the reciprocal of the
rate with the reciprocal of the substrate concentration, the slope will be
$\frac{K_m}{V_m}$ and the intercept will be
$\frac{1}{V_m}$ - giving both $V_m$ and $K_m$.  

The {\em rate} of the reaction is as we have seen a {\em derivative} (slope)
and is $-\frac{dS}{dt}$ (change of substrate concentration with time) or 
$\frac{dP}{dt}$ (change of product with time). 
In the para nitrophenol phosphate (PNPP) to 
para nitrophenol phosphate (PNP) system catalyzed by alkaline
phosphatase, there is a 
one to one stoichiometric relationship between the substrate (PNPP, colorless)
and the product (PNP, yellow in color) and so 
$-\frac{dS}{dt} = \frac{dP}{dt}$.
The substrate concentration at any time $t$ is simply
$S = S0 - P$ where $S0$ is the initial
concentration of the substrate and $P$ the product concentration.
Typically, absorbance values are collected as a function of time (at a wavelength of 
405 nm, where
the product, $P$ (PNP), absorbs) and converted to concentration using the
known value of the extinction coefficient of PNP.
Thus, the data collected in the laboratory is typically $P$ versus time ($t$) 
and so calculation of the {\em rate} requires that a 
{\em derivative} be calculated from the data.  Calculating {\em derivatives} can be
tricky and subject to noise - so one strategy used is to 
calculate the {\em initial} rate (or the derivative at time zero) and that is
far easy to determine.
And thus, a plot of the reciprocal of the initial rate versus the 
reciprocal of the initial substrate concentration can be used to
calculate the $V_m$ and $K_m$ (as shown in equation~\ref{lwbequation}.

\section{Using initial rate data at varying initial substrate concentrations}
Consider a data set that 
had 24 experiments - four different starting enzyme
concentrations and six different starting substrate concentrations.  The concentrations
and the names of the files are given in the table~\ref{tab:dataset}

\begin{table}
\resizebox{\textwidth}{!}{
\begin{tabular}{|c|c|c|c|c|c|c|}
\cline{1-7}
Enzyme Concentration & \multicolumn{6}{ c| }{Substrate Concentration} \\ 
(moles/liter) & \multicolumn{6}{ c| }{moles/liter} \\ \cline{2-7}
 & 0.0003 & 0.00015 & 0.000075 & 0.00003 & 0.000015 & 0.0000075 \\ \cline{1-7}
15.36 x 10$^{-9}$ & zs1e1.dat & zs2e1.dat & zs3e1.dat & zs4e1.dat & zs5e1.dat & zs6e1.dat \\ \cline{1-7}
7.68 x 10$^{-9}$ & zs1e2.dat & zs2e2.dat & zs3e2.dat & zs4e2.dat & zs5e2.dat & zs6e2.dat \\ \cline{1-7}
3.84 x 10$^{-9}$ & zs1e3.dat & zs2e3.dat & zs3e3.dat & zs4e3.dat & zs5e3.dat & zs6e3.dat \\ \cline{1-7}
1.92 x 10$^{-9}$ & zs1e3.dat & zs2e4.dat & zs3e4.dat & zs4e4.dat & zs5e4.dat & zs6e4.dat \\ \cline{1-7}
\hline
\end{tabular}
}
\caption{Data sets at different starting enzyme and substrate concentrations}
\label{tab:dataset}
\end{table}

<<initialrates.py>>=
<<requiredlibraries>>
<<readnonblanklines>>

def readdata(datafile,extinction):
    X = []
    Y1 = []
    Y = []
    myfile = open(datafile,'r')
    for line in nonblank_lines(myfile):
        X.append(float(line.split()[0]))
        Y1.append(float(line.split()[1]))
    for i in range(len(Y1)):
        Y.append((Y1[i] - Y1[0])/extinction)
    return X, Y

def readlistoffilenames(listoffiles):
    filenames = []
    S = []
    myfile = open(listoffiles,'r')
    for line in nonblank_lines(myfile):
        filenames.append(str(line.split()[0]))
        S.append(float(line.split()[1]))
    return filenames, S    

def lineweaverburke(listoffiles,extinction,E0):
    myfiles = readlistoffilenames(listoffiles)[0]
    mysubs = readlistoffilenames(listoffiles)[1]
    lwbx = []
    lwby = []
    for i in range(len(myfiles)):
        filename = myfiles[i]
        sconc = mysubs[i]
        X = readdata(filename,extinction)[0]
        Y = readdata(filename,extinction)[1]
        polycoeffs = scipy.polyfit(X[0:3],Y[0:3],1) 
        totalpolyfitcoeffs = scipy.polyfit(X,Y,3)
        totalyfit = scipy.polyval(totalpolyfitcoeffs,X)
        yfit = scipy.polyval(polycoeffs,X)  
        plt.plot(X,Y)
        plt.plot(X,yfit)
        plt.plot(X,totalyfit)
#        plt.show()
        plt.clf()
        mypoly = np.poly1d(np.array(totalpolyfitcoeffs))
        myderiv = np.polyder(mypoly)
        lwbx.append(1.0/sconc)
        lwby.append(1.0/myderiv(X[0]))

    
    lwbpolycoeffs = scipy.polyfit(lwbx,lwby,1) 
    lwbyfit = scipy.polyval(lwbpolycoeffs,lwbx)  
 
    plt.xlabel('Reciprocal of Substrate Concentration')
    plt.ylabel('Reciprocal of Initial Rate of Reaction')
    plt.title('Lineweaver Burke Plot and Linear Fit')
    plt.grid()
    plt.plot(lwbx,lwby)
    plt.plot(lwbx,lwbyfit,linestyle='dashed')
    plt.savefig('lvbplot.png')
    plt.clf()
    slope, intercept, r, p, err = stats.linregress(lwbx,lwby)
    Vm = 1.0/intercept 
    Km = slope/intercept
    k2 = Vm/E0
    Vminit = Vm
    Kminit = Km
    k2init = k2
    fname = 'resultsfile'
    rfile = open(fname, 'a+')
    rfile.write(listoffiles)
    rfile.write('\n')
    rfile.write('Vm = ' + str(Vm) + '\n')
    rfile.write('Km = ' + str(Km) + '\n')
    rfile.write('k2 = ' + str(k2) + '\n')
    rfile.close()
    return Vminit, Kminit, k2init

Vminitrates = []
Kminitrates = []
k2initrates = []
E0initrates = []

listoffiles = 'zconcentration1.dat'
E0 = 15.36e-9
extinction = 18800.0
(Vminit, Kminit, k2init) = lineweaverburke(listoffiles,extinction,E0)

Vminitrates.append(Vminit)
Kminitrates.append(Kminit)
E0initrates.append(E0)
k2initrates.append(k2init)

listoffiles = 'zconcentration2.dat'
E0 = 7.68e-9
extinction = 18800.0
(Vminit, Kminit, k2init) = lineweaverburke(listoffiles,extinction,E0)

Vminitrates.append(Vminit)
Kminitrates.append(Kminit)
E0initrates.append(E0)
k2initrates.append(k2init)

listoffiles = 'zconcentration3.dat'
E0 = 3.84e-9
extinction = 18800.0
(Vminit, Kminit, k2init) = lineweaverburke(listoffiles,extinction,E0)

Vminitrates.append(Vminit)
Kminitrates.append(Kminit)
E0initrates.append(E0)
k2initrates.append(k2init)

listoffiles = 'zconcentration4.dat'
E0 = 1.92e-9
extinction = 18800.0
(Vminit, Kminit, k2init) = lineweaverburke(listoffiles,extinction,E0)

Vminitrates.append(Vminit)
Kminitrates.append(Kminit)
E0initrates.append(E0)
k2initrates.append(k2init)

@

\begin{pycode}
from initialrates import *
\end{pycode}

\begin{table}[!htp]
\begin{tabular}{|c|c|c|c|}
\hline
$V_m$ & $k_2$ & $E_0$ & $K_m$ \\
\hline
\py{Vminitrates[0]}& \py{k2initrates[0]} & \py{E0initrates[0]} & \py{Kminitrates[0]} \\
\py{Vminitrates[1]}& \py{k2initrates[1]} & \py{E0initrates[1]} & \py{Kminitrates[1]} \\
\py{Vminitrates[2]}& \py{k2initrates[2]} & \py{E0initrates[2]} & \py{Kminitrates[2]} \\
\py{Vminitrates[3]}& \py{k2initrates[3]} & \py{E0initrates[3]} & \py{Kminitrates[3]} \\
\hline
\end{tabular}
\end{table}

\subsection{Calculating $K_m$ and $V_m$: Using the Integrated Rate Equation}
\label{one}
The Michaelis Menten rate equation~\ref{michaelismentenrate} can be integrated to derive a
relationship between substrate (or product) concentration and time. 
The equation ~\ref{michaelismentenrate} can be rearranged 

\begin{equation}
-\int_{S0}^{S}\left[\frac{[S]+K_m}{V_m[S]}\right] dS = \int_{0.0}^{t} dt
\label{mmintegrated}
\end{equation}

and 
integrated from initial time $t=0$ when $S = S0$ to time $t$ 
when substrate concentration is $S$
The integrated form of the Michaelis Menten rate equation is 

\begin{equation}
\frac{S - S0}{V_m} + \frac{K_m}{V_m}\mbox{ln}\left(\frac{S}{S0}\right) = - t
\label{mmintegrated1}
\end{equation}

and changing signs, we get 

\begin{equation}
\frac{S0 - S}{V_m} + \frac{K_m}{V_m}\mbox{ln}\left(\frac{S0}{S}\right) = t
\label{integratedform1}
\end{equation}

The integrated form, equation~\ref{integratedform1}~can be used to calculate the $S$ at any time $t$ given the
initial substrate concentration (knowing $V_m$ and $K_m$)
OR to calculate the time given the 
substrate concentration at any time and its initial concentration.

On the other hand, if we had experimental data on 
$S$ verus $t$, we can calculate the {\em best fit} values
for $K_m$ and $V_m$. 
Instead of expressing $t$ as a function of $S$ (substrate concentration) we can convert
it in terms of $P$ (product) (because we know that there is a 1:1 stoichiometric relationship
between the substrate and product (as shown in equation~\ref{mechanism1}). 

We rearrange equation~\ref{integratedform1} by noting that 
$S0 - S = P$ and so $S = S0 - P$, we get

\begin{equation}
\frac{P}{V_m} + \frac{K_m}{V_m}\mbox{ln}\left(\frac{S0}{S0 - P}\right) = t
\label{mmintegrated}
\end{equation}

Experimental data that has $P$ versus time ($t$) can be used to determine the
best values of $K_m$ and $V_m$ (without having to do the experiments at different
starting substrate concentrations, to get the initial rates).
The initial substrate concentration is assumed to be known for each data set.
The product concentration ($P$) is calculated from the
absorbance values using a known value of the extinction coefficient,
reported as 18,800.0 for the para nitrophenol.  Since we do not know the
precise value of this extinction coefficient, in the analysis of the data sets
using the integral forms as in ~\ref{integratedform1} 
we can if we choose to, allow the $S0$ to vary
and obtaine values of $V_m$ and $K_m$.

\subsection{Allowing $S0$, $V_m$ and $K_m$ to vary}

<<directfitPvsTime.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
def sintegrated(x,Km, Vm, S0):
    return (Km/Vm)*np.log(S0/(S0 - x)) + x/Vm
infile = 'zs6e4.dat'
(Td, Zd) = abs2conc(infile)

popt, pcov = curve_fit(sintegrated, Zd, Td, p0 = (0.001, 0.001, 0.00075))

#print popt

myKm1 = popt[0]
myVm1 = popt[1]
myS01 = popt[2]

myTd = []
for i in range(len(Zd)):
    myTd.append((myKm1/myVm1)*np.log(myS01/(myS01 - Zd[i])) + Zd[i]/myVm1)
    
plt.plot(Td,Zd)
plt.plot(myTd, Zd)
plt.savefig('directfit.png')
plt.clf()
#plt.show()

@



\subsection{Allowing $V_m$ and $K_m$ to vary, $S0$ fixed}


<<directfitPvsTimeSfixed.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
def sintegrated(x,Km, Vm):
    return (Km/Vm)*np.log(S0/(S0 - x)) + x/Vm

infile = 'zs6e4.dat'
(Td, Zd) = abs2conc(infile)

S0 = 0.0000075
S02 = S0
#S0 = 6.136e-6

popt, pcov = curve_fit(sintegrated, Zd, Td, p0 = (0.001, 0.001))

#print popt

myKm2 = popt[0]
myVm2 = popt[1]

myTd = []
for i in range(len(Zd)):
    myTd.append((myKm2/myVm2)*np.log(S02/(S02 - Zd[i])) + Zd[i]/myVm2)
    
plt.plot(Td,Zd)
plt.plot(myTd, Zd)
plt.savefig('directfitS0fixed.png')
plt.clf()
#plt.show()

@

\begin{pycode}
from directfitPvsTime import *
\end{pycode}

\begin{table}[!htp]
\begin{tabular}{|c|c|c|}
\hline
$V_m$ & $K_m$ & $S_0$ \\
\hline
\py{myVm1} & \py{myKm1} & \py{myS01} \\
\hline
\end{tabular}
\caption{S0 allowed to vary}
\end{table}

\begin{pycode}
from directfitPvsTimeSfixed import *
\end{pycode}

\begin{table}[!htp]
\begin{tabular}{|c|c|c|}
\hline
$V_m$ & $K_m$ & $S_0$ \\
\hline
\py{myVm2} & \py{myKm2} & \py{S02}  \\
\hline
\end{tabular}
\caption{S0 Fixed}
\end{table}

%\py{Vmc}, \py{Kmc}, \py{Vmi}, \py{Kmi}, 
%\py{Vmcg}, \py{nVmc}, \py{Kmcg}, \py{nKmc}, \py{Vmig}, 
%\py{nVmi}, \py{Kmig}, \py{nKmi}

\section{The Enzyme Mechanism, complete}
A mechanism proposed for the conversion of a substrate 
$PNPP$ (para nitrophenol phosphate)
to a product, $PNP$ (para nitrophenol) is shown in ~\ref{mechanism1}.  The substrate is postulated to bind
to the enzyme (E) forming an enzyme 
substrate complex (E--PNPP) which then breaks up into the original enzyme
and the product (PNP), with the release of free phosphate (P).
This very familiar reaction can be written like this (where we show ALL of the
participants in the reaction)

\begin{equation}
\ce{E + PNPP <=>T[\cf{k_1}][\cf{k_{-1}}] E--PNPP ->T[\cf{k_{2}}] E + PNP + P}
\label{mechanism1}
\end{equation}
This reaction is often written using $S$ for substrate
(shown as PNPP in ~\ref{mechanism1} 
and $P$ for the product (shown as PNP in ~\ref{mechanism1}).
We will use $Ph$ for the free phosphate 

\begin{equation}
\ce{E + S <=>T[\cf{k_1}][\cf{k_{-1}}] ES ->T[\cf{k_{2}}] E + P + Ph}
\label{mechanism2}
\end{equation}

The free phosphate shown in ~\ref{mechanism2} can bind to the
enzyme as a competitive inhibitor, so we can imagine one
additional reaction

\begin{equation}
\ce{E + Ph <=>T[\cf{k_3}][\cf{k_{-3}}] E(Ph)}
\label{mechanism3}
\end{equation}

In a batch reactor, we can write expressions for the rate
of change for the different species involved in the reactions
shown in ~\ref{mechanism2} and in \ref{mechanism3}.

\begin{center}
\begin{subequations}\label{completeset}
\begin{align}
\frac{d[S]}{dt} & = & -k_1[E][S] + k_{-1}[ES] \label{completeset:S} \\
\frac{d[E]}{dt} & = & -k_1[E][S] + k_{-1}[ES] + k_2[ES] - k_{3}[E][Ph] + k_{-3}[EPh] \label{completeset:E} \\
\frac{d[ES]}{dt} & = & k_1[E][S]- k_{-1}[ES] - k_2[ES] \label{completeset:ES} \\
\frac{d[P]}{dt} & = & k_2[ES] \label{completeset:P} \\
\frac{d[EPh]}{dt} & = & k_3[E][Ph] - k_{-3}[EPh] \label{completeset:EPh} \\
\frac{d[Ph]}{dt} & = & k_{-3}[EPh] - k_{3}[E][Ph] \label{completeset:Ph} 
\end{align}
\end{subequations}
\end{center}

This system of differential equations shown as
~\ref{completeset} can be integrated if we know the rate constants
and initial values of the species.  The python code to do this is shown
as a function called {\em completeset} (and we will see how to use this
function later)

<<completeset>>=

def completeset(parameters):

    k1, km1, k2, k3, km3, S0, E0, ES0, P0, EPh0, Ph0, start_time, end_time, intervals = parameters
    rates = (k1, km1, k2, k3, km3)
    initial_cond = [S0, E0, ES0, P0, EPh0, Ph0]
# define the differential equation    
    def eq(par,initial_cond,start_t,end_t,incr):
#-time-grid-----------------------------------
      t  = np.linspace(start_t, end_t,incr)
#differential-eq-system----------------------
      def funct(y,t):
        S=y[0]
        E=y[1]
        ES=y[2]
        P=y[3]
        EPh=y[4]
        Ph=y[5]
        k1,km1,k2,k3,km3=par
        # the model equations (see Munz et al. 2009)
        f0 = -k1*E*S + km1*ES 
        f1 = -k1*E*S + km1*ES + k2*ES -k3*E*Ph + km3*EPh
        f2 = k1*E*S - km1*ES - k2*ES
        f3 = k2*ES
        f4 = k3*E*Ph - km3*EPh
        f5 = -k3*E*Ph + km3*EPh
        return [f0, f1, f2, f3, f4, f5]
#integrate------------------------------------
      ds = integrate.odeint(funct,initial_cond,t)
      return (ds[:,0],ds[:,1],ds[:,2],ds[:,3],ds[:,4],ds[:,5],t)

# initial conditions for the calculations
    y0 = [S0, E0, ES0, P0, EPh0, Ph0]
# model steps
#---------------------------------------------------

# F0 = S, F1 =E, F2 = ES, F3 = P, F4 = EPh, F5 = Ph
    F0,F1,F2,F3,F4,F5,T=eq(rates,y0,start_time,end_time,intervals)
# calculate Vm using the k2 and E0 
    Vm = k2*E0
# calculate the Km using k1, k2 and km1    
    Km = (km1 + k2)/k1
# you can calculate the equilibrium constant at all the time points    
    Keq = []
    for i in range(len(T)):
        Keq.append(F2[i]/(F0[i]*F1[i]))
    return F0, F1, F2, F3, F4, F5, T, Keq

@


\section{The Enzyme Mechanism, approximate}

If we ignore the competitive binding
of phosphate to the enzyme, we get a simpler set of equations shown in
~\ref{approximateset} (we can also get this simpler set by 
setting $k_3$ and $k_{-3}$ to zero)

\begin{center}
\begin{subequations}\label{approximateset}
\label{diffeqs}
\begin{align}
\frac{d[S]}{dt} & = & -k_1[E][S] + k_{-1}[ES] \label{approximateset:S} \\
\frac{d[E]}{dt} & = & -k_1[E][S] + k_{-1}[ES] + k_2[ES] \label{approximateset:E} \\
\frac{d[ES]}{dt} & = & k_1[E][S]- k_{-1}[ES] - k_2[ES] \label{approximateset:ES} \\
\frac{d[P]}{dt} & = & k_2[ES] \label{approximateset:P} 
\end{align}
\end{subequations}
\end{center}

This system of differential equations shown as
~\ref{approximateset} can also be integrated if we know the rate constants
and initial values of the species.  The python code to do this is shown
as a function called {\em approximateset} (and we will see how to use this
function later)


<<approximateset>>=
def approximateset(parameters):

    k1, km1, k2, S0, E0, ES0, P0, start_time, end_time, intervals = parameters
    rates = (k1, km1, k2)
    initial_cond = [S0, E0, ES0, P0]
# define the differential equation    
    def eq(par,initial_cond,start_t,end_t,incr):
#-time-grid-----------------------------------
      t  = np.linspace(start_t, end_t,incr)
#differential-eq-system----------------------
      def funct(y,t):
        S=y[0]
        E=y[1]
        ES=y[2]
        P=y[3]
        k1,km1,k2=par
        # the model equations (see Munz et al. 2009)
        f0 = -k1*E*S + km1*ES 
        f1 = -k1*E*S + km1*ES + k2*ES
        f2 = k1*E*S - km1*ES - k2*ES
        f3 = k2*ES
        return [f0, f1, f2, f3]
#integrate------------------------------------
      ds = integrate.odeint(funct,initial_cond,t)
      return (ds[:,0],ds[:,1],ds[:,2],ds[:,3],t)

# initial conditions for the calculations
    y0 = [S0, E0, ES0, P0]
# model steps
#---------------------------------------------------

# F0 = S, F1 =E, F2 = ES, F3 = P
    F0,F1,F2,F3,T=eq(rates,y0,start_time,end_time,intervals)
# calculate Vm using the k2 and E0 
    Vm = k2*E0
# calculate the Km using k1, k2 and km1    
    Km = (km1 + k2)/k1
# you can calculate the equilibrium constant at all the time points    
    Keq = []
    for i in range(len(T)):
        Keq.append(F2[i]/(F0[i]*F1[i]))
    return F0, F1, F2, F3, T, Keq

@

\section{Fit data to complete set}
If we have data for example on how the product changes with time, we can
{\em fit} the data to {\em estimate} the rate constants.
The python function to fit the completeset 
of equations is called fitcompleteset.

<<fitcompleteset>>=
def fitcompleteset(Td, Zd, parameters):

#Td, Zd contain - time and concentration (time is in seconds)

    k1, km1, k2, k3, km3, S0, E0, ES0, P0, EPh0, Ph0, start_time, end_time, intervals = parameters
    rates = (k1, km1, k2, k3, km3)
# define the differential equation    
    def eq(par,initial_cond,start_t,end_t,incr):
#-time-grid-----------------------------------
        t  = np.linspace(start_t, end_t,incr)
#differential-eq-system----------------------
        def funct(y,t):
            S=y[0]
            E=y[1]
            ES=y[2]
            P=y[3]
            EPh=y[4]
            Ph=y[5]
            k1,km1,k2,k3,km3 = par         
        # the model equations f0 is dS/dt, f1 is dE/dt, f2 is dES/dt, F3 is dP/dt
            f0 = -k1*E*S + km1*ES 
            f1 = -k1*E*S + km1*ES + k2*ES
            f2 = k1*E*S - km1*ES - k2*ES
            f3 = k2*ES
            f4 = k3*E*Ph - km3*EPh
            f5 = -k3*E*Ph + km3*EPh
            return [f0, f1, f2, f3, f4, f5]

#integrate - solution returned in order for S, E, ES and P 
        ds = integrate.odeint(funct,initial_cond,t)
        return (ds[:,0],ds[:,1],ds[:,2],ds[:,3],ds[:,4], ds[:,5], t)
# initial conditions for the calculations
    y0 = [S0, E0, ES0, P0, EPh0, Ph0]
# model steps
#---------------------------------------------------
    mt=np.linspace(start_time,end_time,intervals)
 
# model index to compare to data
#----------------------------------------------------
    findindex=lambda x:np.where(mt>=x)[0][0]
    mindex=map(findindex,Td)
#=======================================================
 
#Score Fit of System
#=========================================================
    def score(parms):
#a.Get Solution to system
     S, E, ES, P, EPh, Ph, T=eq(parms,y0,start_time,end_time,intervals)
#b.Pick of Model Points to Compare - Note we pick the product P
     Zm=P[mindex]
#c.Score Difference between model and data points
     ss=lambda data,model:((data-model)**2).sum()
     return ss(Zd,Zm)
#========================================================
 
 
#Optimize Fit
#=======================================================
    fit_score=score(rates)
    answ=fmin(score,(rates),full_output=1,maxiter=100000,maxfun=100000)
    bestrates=answ[0]
    bestscore=answ[1]
    nk1, nkm1, nk2, nk3, nkm3 =answ[0]
    newrates=(nk1, nkm1, nk2, nk3, nkm3)
#=======================================================
 
#Generate optimized Solution to System
#=======================================================
    S, E, ES, P, EPh, Ph, T=eq(newrates,y0,start_time,end_time,intervals)
    Zm=P[mindex]
    Tm=T[mindex]

    Vm = nk2*E0
    Km = (nkm1 + nk2)/nk1
    return T,P,Tm,Zm,nk1,nkm1,nk2, nk3, nkm3

@

\section{Fit data to approximate set}
The python function to fit the approximate set of 
differential equations is called fitapproximateset

<<fitapproximateset>>=
def fitapproximateset(Td, Zd, parameters):

#Td, Zd contain - time and concentration (time is in seconds)

    k1, km1, k2, S0, E0, ES0, P0, start_time, end_time, intervals = parameters
    rates = (k1, km1, k2)
# define the differential equation    
    def eq(par,initial_cond,start_t,end_t,incr):
#-time-grid-----------------------------------
        t  = np.linspace(start_t, end_t,incr)
#differential-eq-system----------------------
        def funct(y,t):
            S=y[0]
            E=y[1]
            ES=y[2]
            P=y[3]
            k1,km1,k2=par
        # the model equations f0 is dS/dt, f1 is dE/dt, f2 is dES/dt, F3 is dP/dt
            f0 = -k1*E*S + km1*ES 
            f1 = -k1*E*S + km1*ES + k2*ES
            f2 = k1*E*S - km1*ES - k2*ES
            f3 = k2*ES
            return [f0, f1, f2, f3]
#integrate - solution returned in order for S, E, ES and P 
        ds = integrate.odeint(funct,initial_cond,t)
        return (ds[:,0],ds[:,1],ds[:,2],ds[:,3],t)
# initial conditions for the calculations
    y0 = [S0, E0, ES0, P0]
# model steps
#---------------------------------------------------
    mt=np.linspace(start_time,end_time,intervals)
 
# model index to compare to data
#----------------------------------------------------
    findindex=lambda x:np.where(mt>=x)[0][0]
    mindex=map(findindex,Td)
#=======================================================
 
#Score Fit of System
#=========================================================
    def score(parms):
#a.Get Solution to system
     S, E, ES, P, T=eq(parms,y0,start_time,end_time,intervals)
#b.Pick of Model Points to Compare - Note we pick the product P
     Zm=P[mindex]
#c.Score Difference between model and data points
     ss=lambda data,model:((data-model)**2).sum()
     return ss(Zd,Zm)
#========================================================
 
 
#Optimize Fit
#=======================================================
    fit_score=score(rates)
    answ=fmin(score,(rates),full_output=1,maxiter=100000,maxfun=100000)
    bestrates=answ[0]
    bestscore=answ[1]
    nk1, nkm1, nk2 =answ[0]
    newrates=(nk1, nkm1, nk2)
#=======================================================
 
#Generate optimized Solution to System
#=======================================================
    S, E, ES, P, T=eq(newrates,y0,start_time,end_time,intervals)
    Zm=P[mindex]
    Tm=T[mindex]

    Vm = nk2*E0
    Km = (nkm1 + nk2)/nk1
    return T,P,Tm,Zm,nk1,nkm1,nk2 

@

<<setupsystem>>=
# set the values for rate constants
k1 = 2.912e+09
km1 = 4.275e+05
k2 = 300.0
k3 = 1.835e+08
km3 = 600.0
# set initial conditions
S0 = 7.5e-06
E0 = 1.92e-09
ES0 = 0.0
P0 = 0.0
EPh0 = 0.0
Ph0 = 7.5e-06
Vm = k2*E0
Km = (km1 + k2)/k1 
# set conditions for integration, time in seconds
start_time = 0.0
end_time = 300.0
intervals = 1000

parameters1 = (k1, km1, k2, S0, E0, ES0, P0, start_time, end_time, intervals)
parameters2 = (k1, km1, k2, k3, km3, S0, E0, ES0, P0, EPh0, Ph0, start_time, end_time, intervals)

infile = 'zs6e4.dat'

@


<<completemodel.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
<<completeset>>
<<setupsystem>>
(myS, myE, myES, myproduct, myEPh, myPh, mytime, myKeq) = completeset(parameters2)

plt.xlabel(' Time in Seconds ')
plt.ylabel(' Concentration, moles/liter ')
plt.plot(mytime, myproduct)
plt.title('Complete Model, including effect of phosphate')
plt.savefig('completemodel.png')
plt.clf()
Vmc = Vm
Kmc = Km
#plt.show()

@

<<bothmodels.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
<<completeset>>
<<setupsystem>>
(mySc, myEc, myESc, myproductc, myEPhc, myPhc, mytimec, myKeqc) = completeset(parameters2)
<<approximateset>>
(mySi, myEi, myESi, myproducti, mytimei, myKeqi) = approximateset(parameters1)

plt.xlabel(' Time in Seconds ')
plt.ylabel(' Concentration, moles/liter ')
plt.plot(mytimec, myproductc)
plt.plot(mytimei, myproducti)
plt.title('Both Models')
plt.savefig('bothmodels.png')
plt.clf()
#plt.show()

@

\begin{pycode}
from completemodel import *
\end{pycode}

We run {\em completemodel.py} using the values for the
different rate constants and initial concentrations as shown in the
table

\begin{table}
\begin{tabular}{|c|c|c|c|c|c|}
\hline
$k_1$ & $k_{-1}$ & $k_2$ & $k_3$ & $k_{-3}$ & \\
\hline
\py{k1} & \py{km1} & \py{k2} & \py{k3} & \py{km3} & \\
\hline
$S_0$ & $E_0$ & $ES_0$ & $EPh_0$ & $Ph_0$ & $P_0$ \\
\hline
\py{S0} & \py{E0} & \py{ES0} & \py{EPh0} & \py{Ph0} & \py{P0} \\
\hline
$V_m$ & $K_m$ & & & & \\
\hline
\py{Vmc} & \py{Kmc} & & & & \\
\hline
 
\end{tabular}
\caption{The Complete Model}
\end{table}

<<approximatemodel.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
<<approximateset>>
<<setupsystem>>
(myS, myE, myES, myproduct, mytime, myKeq) = approximateset(parameters1)
plt.xlabel(' Time in Seconds ')
plt.ylabel(' Concentration, moles/liter ')
plt.plot(mytime, myproduct)
plt.title('Approximate Model, no phosphate inhibition')
plt.savefig('approximatemodel.png')
plt.clf()
Vmi = Vm
Kmi = Km

@




\begin{pycode}
from approximatemodel import *
\end{pycode}

We run {\em approximatemodel.py} using the values for the
different rate constants and initial concentrations as shown in the
table

\begin{table}
\begin{tabular}{|c|c|c|c|c|c|}
\hline
$k_1$ & $k_{-1}$ & $k_2$ & & &  \\
\hline
\py{k1} & \py{km1} & \py{k2} & & & \\
\hline
$S_0$ & $E_0$ & $ES_0$ & & & $P_0$ \\
\hline
\py{S0} & \py{E0} & \py{ES0} & & & \py{P0} \\
\hline
$V_m$ & $K_m$ & & & & \\
\hline
\py{Vmi} & \py{Kmi} & & & & \\
\hline
\caption{The approximate model}
\end{tabular}
\end{table}


<<fitcompletemodel.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
<<approximateset>>
<<setupsystem>>
<<completeset>>
<<fitcompleteset>>

(Td, Zd) = abs2conc(infile)

(Tfit,Pfit,Tmfit,Zmfit,nk1,nkm1,nk2, nk3, nkm3) = fitcompleteset(Td,Zd,parameters2)

plt.plot(Tfit,Pfit)

parameters2 = (nk1, nkm1, nk2, nk3, nkm3, S0, E0, ES0, P0, EPh0, Ph0, start_time, end_time, intervals)

(mySf, myEf, myESf, myproductf, myEPh, myPh, mytimef, myKeq) = completeset(parameters2)

Vmcg = k2*E0
Kmcg = (km1 + k2)/k1

nVmc = nk2*E0
nKmc = (nkm1 + nk2)/nk1 

plt.plot(mytimef,myproductf,'ro',markevery=40)
#plt.show()
plt.savefig('fitcompletemodel.png')
plt.clf()

@



\begin{pycode}
from fitcompletemodel import *
\end{pycode}

\begin{table}
\begin{tabular}{|c|c|c|c|c|c|}
\hline
$k_1$ & $k_{-1}$ & $k_2$ & $k_3$ & $k_{-3}$ & \\
\hline
\py{k1} & \py{km1} & \py{k2} & \py{k3} & \py{km3} & \\
\hline
$k_1$ & $k_{-1}$ & $k_2$ & $k_3$ & $k_{-3}$ & \\
\hline
\py{nk1} & \py{nkm1} & \py{nk2} & \py{nk3} & \py{nkm3} & \\
\hline
\py{S0} & \py{E0} & \py{ES0} & \py{EPh0} & \py{Ph0} & \py{P0} \\
\hline
$V_m$ (guess) & $V_m$ (final) & $K_m$ (guess) & $K_m$ (final) & & \\
\hline
\py{Vmcg} & \py{nVmc} & \py{Kmcg} & \py{nKmc} & & \\
\hline 
\end{tabular}
\caption{Fitting the complete Model}
\end{table}


<<fitapproximatemodel.py>>=
<<requiredlibraries>>
<<readnonblanklines>>
<<absorbance2concentration>>
<<approximateset>>
<<setupsystem>>
<<approximateset>>
<<fitapproximateset>>

(Td, Zd) = abs2conc(infile)

(Tfit,Pfit,Tmfit,Zmfit,nk1,nkm1,nk2) = fitapproximateset(Td,Zd,parameters1)

plt.plot(Tfit,Pfit)

parameters1 = (nk1, nkm1, nk2, S0, E0, ES0, P0, start_time, end_time, intervals)

(mySf, myEf, myESf, myproductf, mytimef, myKeq) = approximateset(parameters1)

Vmig = k2*E0
Kmig = (km1 + k2)/k1

nVmi = nk2*E0
nKmi = (nkm1 + nk2)/nk1 


plt.plot(mytimef,myproductf,'ro',markevery=40)
plt.savefig('fitapproximatemodel.png')
plt.clf()
#plt.show()



@




\begin{pycode}
from fitapproximatemodel import *
\end{pycode}

\begin{table}
\begin{tabular}{|c|c|c|c|c|c|}
\hline
$k_1$ & $k_{-1}$ & $k_2$ & & & \\
\hline
\py{k1} & \py{km1} & \py{k2} & &  & \\
\hline
$k_1$ & $k_{-1}$ & $k_2$ & & & \\
\hline
\py{nk1} & \py{nkm1} & \py{nk2} &  & & \\
\hline
\py{S0} & \py{E0} & \py{ES0} & \py{EPh0} & \py{Ph0} & \py{P0} \\
\hline
$V_m$ (guess) & $V_m$ (final) & $K_m$ (guess) & $K_m$ (final) & & \\
\hline
\py{Vmig} & \py{nVmi} & \py{Kmig} & \py{nKmi} & & \\
\hline
\end{tabular}
\caption{Fitting the approximate model}
\end{table}


%\lstinputlisting{completemodel.py}
%\section{python listing for fitcompletemodel}
%\lstinputlisting[language=python]{fitcompletemodel.py}
%\section{python listing for fitapproximatemodel}
%\lstinputlisting[language=python]{fitapproximatemodel.py}



<<studyenzymekinetics.make>>=
all: study-enzyme-kinetics.pdf 

study-enzyme-kinetics.tex: study-enzyme-kinetics.w
	noweb study-enzyme-kinetics.w

completemodel.png: completemodel.py
	python2 completemodel.py

approximatemodel.png: approximatemodel.py
	python2 approximatemodel.py

fitapproximatemodel.png: fitapproximatemodel.py
	python2 fitapproximatemodel.py

directfit.png: directfitPvsTime.py
	python2 directfitPvsTime.py

directfitS0fixed.png: directfitPvsTimeSfixed.py
	python2 directfitPvsTimeSfixed.py

bothmodels.png: bothmodels.py
	python2 bothmodels.py

lvbplot.png: initialrates.py
	python2 initialrates.py

study-enzyme-kinetics.pdf: study-enzyme-kinetics.tex completemodel.png approximatemodel.png fitapproximatemodel.png directfit.png directfitS0fixed.png bothmodels.png lvbplot.png

	pdflatex -shell-escape study-enzyme-kinetics.tex
	pythontex --interpreter python:python2 study-enzyme-kinetics.tex
	pdflatex -shell-escape study-enzyme-kinetics.tex
	pdflatex -shell-escape study-enzyme-kinetics.tex

clean:
	rm -r -f *.aux *.toc study-enzyme-kinetics.pdf study-enzyme-kinetics.tex _minted*
	rm -r -f pythontex-files-study-enzyme-kinetics
	rm -f study-enzyme-kinetics.pytxcode

@


\end{document}

