% I have used noweb as the literate programming tool
% noweb nonlinear-regression-example.nw will extract the files you need
% the files include: a make file, a python file, a tex file 
% You will need a tex distribution (on windows, you can use www.miktex.org)
% and python (anaconda will work nicely) and pythontex (to help include python within tex)
% the make file is called linreg.make
% and on unix systems, if you type make -f linreg.make you will get the pdf 

\documentclass{tufte-handout}
\usepackage{noweb}
\usepackage{hyperref}
\usepackage{graphicx}
\usepackage{amsmath}  % extended mathematics
\usepackage{pythontex}
\setpythontexworkingdir{.}
\hypersetup
{   pdfauthor = {Krishnan Chittur},
  pdftitle={An example in nonlinear regression},
  colorlinks=TRUE,
  linkcolor=black,
  citecolor=blue,
  urlcolor=blue
}

\title{Nonlinear Regression}
\author{Krishnan Chittur \\ \url{http://webpages.uah.edu/~chitturk}}
\date{January, 2018}

\begin{document}
\tableofcontents

\section{Abstract}
In this document, I show how you can use
{\em literal programming} - where you can 
mix documentation, code and sample data files together.
A program called
{\em noweb} can be used to create/extract the documentation and the code.
You can embed commands to execute the code and include results from running the
code in the document.
You will need \href{https://www.cs.tufts.edu/~nr/noweb/}{noweb} if you begin with this
file. The command {\em noweb nonlinear-regression-example.w} will extract the
\TeX file and the python codes.  You can then run
{\em make -f linreg.make} to {\bf make} the documentation and run the code.

<<nonlinreg.make>>=
all: nonlinear-regression-example.pdf 

nonlinear-regression-example.tex: nonlinear-regression-example.w
	noweb nonlinear-regression-example.w
figure.png: dononlinear.py
	python2 dononlinear.py
noisyfigure.png: dononlinear.py
	python2 dononlinear.py
nonlinear-regression-example.pdf: nonlinear-regression-example.tex figure.png noisyfigure.png
	pdflatex -shell-escape nonlinear-regression-example.tex
	pythontex --interpreter python:python2 nonlinear-regression-example.tex
	pdflatex -shell-escape nonlinear-regression-example.tex
	pdflatex -shell-escape nonlinear-regression-example.tex

clean:
	rm -r -f figure.png noisyfigure.png *.aux *.toc nonlinear-regression-example.pdf nonlinear-regression-example.tex _minted*
	rm -r -f pythontex-files-nonlinear-regression-example
	rm -f nonlinear-regression-example.pytxcode

@

\subsection{Required Libraries}
A collection of python libraries that many of the scripts will need 

<<requiredlibraries>>=
import numpy as np
import math
import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl import load_workbook
from random import uniform as randunif
import numpy as np
from scipy.optimize import curve_fit

@

\section{A simple nonlinear model}
Imagine you want to fit a data set to a nonlinear model.

\begin{equation}
y = \frac{Vm x}{Km +x}
\end{equation}

I will first simply {\em create} a data set using some value
for the Vm and the Km, add some random noise to the
$y$ value and create an excel file.  

<<createdata>>=
def createdata(Vm, Km, xmin, xmax, addnoise):
    myx = np.linspace(xmin, xmax, 50)
    myy = []
    for i in range(len(myx)):
        if (addnoise == 1):
            myrandom = np.random.normal(-0.1*myx[i],0.1*myx[i])
        else:
            myrandom = 0.0
        myy.append(Vm*(myx[i]+myrandom)/(Km + myx[i]+myrandom))
    return myx, myy

@

I will then create an xlsx file where I will store some data 
using the createdata function

<<createxlsxfile>>=
# this function will create an xlsx file using the createdata function
def createxlsx(Vm, Km, xmin, xmax, addnoise, outputfilename):
    myx, myy = createdata(Vm, Km,xmin, xmax, addnoise)
    outputfile = Workbook()
    mytitle = 'Data'
    mysheet = outputfile.create_sheet(title = mytitle)
    mydetails = outputfile.create_sheet(title = 'Details of Experiment')
    mydetails.cell(row=1,column=1).value = 'Vm'
    mydetails.cell(row=1,column=2).value = 'Km'
    mydetails.cell(row=2,column=1).value = Vm
    mydetails.cell(row=2,column=2).value = Km
    mysheet.cell(row=1,column=1).value = 'X'
    mysheet.cell(row=1,column=2).value = 'Y'

    for i in range(len(myy)):
        mysheet.cell(row=i+2,column=1).value = myx[i]
        mysheet.cell(row=i+2,column=2).value = myy[i]

    std=outputfile.get_sheet_by_name('Sheet')
    outputfile.remove_sheet(std)
    outputfile.save(outputfilename)

@

The readxlsxfile function will read data from an xlsx file 

<<readxlsxfile>>=
def readxlsx(myfilename):
    wb = load_workbook(filename = myfilename)
    datasheet = wb.get_sheet_by_name('Data')
    nr = datasheet.max_row
    nc = datasheet.max_column
    myx = []
    myy = []
    for i in range(nr-1):
        myx.append(datasheet.cell(row=i+2,column=1).value)
        myy.append(datasheet.cell(row=i+2,column=2).value)
    runsheet = wb.get_sheet_by_name('Details of Experiment')
    slope = runsheet.cell(row=2,column=1).value
    intercept = runsheet.cell(row=2,column=2).value
    return myx, myy, Vm, Km

@

nonlinear is a function that returns the value of the nonlinear function
given slope and intercept and is used by the curve fit function to
calculate the best values of slope and intercept

<<nonlinear>>=
def nonlinear(x, Vm, Km):
    return Vm*x/(Km + x)

@

Assemble the complete dononlinear script from the individual chunks

<<dononlinear.py>>=
<<requiredlibraries>>
<<createdata>>
<<createxlsxfile>>
<<readxlsxfile>>
<<nonlinear>>
Vm = 5.0
Km = 1.0
xmin = 0.1
xmax = 10.0
addnoise = 0.0
createxlsx(Vm, Km, xmin, xmax, addnoise, 'file.xlsx')
addnoise = 1.0
createxlsx(Vm, Km, xmin, xmax, addnoise, 'noisyfile.xlsx')

filename = 'file.xlsx'
myx, myy, Vm, Km = readxlsx(filename)
popt, pcov = curve_fit(nonlinear, myx, myy, p0 = (1.0, 1.0))
myVm = popt[0]
myKm = popt[1]
myfittedy = []
for i in range(len(myx)):
	myfittedy.append(myVm*myx[i]/(myKm + myx[i]))

plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.plot(myx, myy,linestyle=':')
plt.plot(myx, myfittedy)
plt.savefig('figure.png')
plt.clf()

filename = 'noisyfile.xlsx'
myx, myy, Vm, Km = readxlsx(filename)
popt, pcov = curve_fit(nonlinear, myx, myy, p0 = (1.0, 1.0))
myVmnoisy = popt[0]
myKmnoisy = popt[1]
myfittedynoisy = []

for i in range(len(myx)):
	myfittedynoisy.append(myVmnoisy*myx[i]/(myKmnoisy + myx[i]))

plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.plot(myx, myy,linestyle=':')
plt.plot(myx, myfittedy)
plt.savefig('noisyfigure.png')
plt.clf()


@

Now execute the dononlinear script 

\begin{pycode}
from dononlinear import *
\end{pycode}

You can access the results/variables from the dononlinear script

\begin{figure}
\includegraphics[scale=0.75]{figure.png}
\caption{No noise added}
\end{figure}

\begin{figure}
\includegraphics[scale=0.75]{noisyfigure.png}
\caption{Noise added}
\end{figure}


The equation was $y = \frac{Vm x}{Km + x}$ The $Vm$ value used
was \py{Vm} and $Km$ value was \py{Km}.

If noise was NOT added to the $y$, the 
calculated slope was \py{myVm} and the 
calculated intercept was \py{myKm}

Some noise was added to the $y$ and when curve fit was used, the 
calculated slope was \py{myVmnoisy} and the 
calculated intercept was \py{myKmnoisy}

\end{document}

