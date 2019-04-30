% I have used noweb as the literate programming tool
% noweb linear-regression-example.nw will extract the files you need
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
  pdftitle={An example in linear regression},
  colorlinks=TRUE,
  linkcolor=black,
  citecolor=blue,
  urlcolor=blue
}

\title{Linear Regression}
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
file. The command {\em noweb linear-regression-example.w} will extract the
\TeX file and the python codes.  You can then run
{\em make -f linreg.make} to {\bf make} the documentation and run the code.

<<linreg.make>>=
all: linear-regression-example.pdf 

linear-regression-example.tex: linear-regression-example.w
	noweb linear-regression-example.w
figure.png: dolinear.py
	python dolinear.py
noisyfigure.png: dolinear.py
	python dolinear.py
linear-regression-example.pdf: linear-regression-example.tex figure.png noisyfigure.png
	pdflatex -shell-escape linear-regression-example.tex
	pythontex --interpreter python:python2 linear-regression-example.tex
	pdflatex -shell-escape linear-regression-example.tex
	pdflatex -shell-escape linear-regression-example.tex

clean:
	rm -r -f figure.png noisyfigure.png *.aux *.toc linear-regression-example.pdf linear-regression-example.tex _minted*
	rm -r -f pythontex-files-linear-regression-example
	rm -f linear-regression-example.pytxcode

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

\section{A simple linear model}
Imagine you want to fit a data set to a linear model.

\begin{equation}
y = slope*x + intercept
\end{equation}

I will first simply {\em create} a data set using some value
for the slope and the intercept, add some random noise to the
$y$ value and create an excel file.  

<<createdata>>=
def createdata(slope, intercept, xmin, xmax, addnoise):
    myx = np.linspace(xmin, xmax, 50)
    myy = []
    for i in range(len(myx)):
        if (addnoise == 1):
            myrandom = np.random.normal(-0.1*myx[i],0.1*myx[i])
        else:
            myrandom = 0.0
        myy.append(slope*(myx[i]+myrandom) + intercept)
    return myx, myy

@

I will then create an xlsx file where I will store some data 
using the createdata function

<<createxlsxfile>>=
# this function will create an xlsx file using the createdata function
def createxlsx(slope, intercept, xmin, xmax, addnoise, outputfilename):
    myx, myy = createdata(slope, intercept, xmin, xmax, addnoise)
    outputfile = Workbook()
    mytitle = 'Data'
    mysheet = outputfile.create_sheet(title = mytitle)
    mydetails = outputfile.create_sheet(title = 'Details of Experiment')
    mydetails.cell(row=1,column=1).value = 'Slope'
    mydetails.cell(row=1,column=2).value = 'Intercept'
    mydetails.cell(row=2,column=1).value = slope
    mydetails.cell(row=2,column=2).value = intercept
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
    return myx, myy, slope, intercept

@

linear is a function that returns the value of the linear function
given slope and intercept and is used by the curve fit function to
calculate the best values of slope and intercept

<<linear>>=
def linear(x,slope, intercept):
    return slope*x + intercept

@

Assemble the complete dolinear script from the individual chunks

<<dolinear.py>>=
<<requiredlibraries>>
<<createdata>>
<<createxlsxfile>>
<<readxlsxfile>>
<<linear>>
from matplotlib import rc 
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

slope = 5.0
intercept = 1.0
xmin = 0.1
xmax = 10.0
addnoise = 0.0
createxlsx(slope, intercept, xmin, xmax, addnoise, 'file.xlsx')

addnoise = 1.0
createxlsx(slope, intercept, xmin, xmax, addnoise, 'noisyfile.xlsx')

filename = 'file.xlsx'
myx, myy, slope, intercept = readxlsx(filename)
popt, pcov = curve_fit(linear, myx, myy, p0 = (1.0, 1.0))
myslope = popt[0]
myintercept = popt[1]
myfittedy = []
for i in range(len(myx)):
	myfittedy.append(myslope*myx[i]+myintercept)

plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.plot(myx, myy,linestyle=':')
plt.plot(myx, myfittedy)
plt.title(r"Equation is y = %s x + %s " %(slope, intercept))   
plt.savefig('figure.png')
plt.clf()


filename = 'noisyfile.xlsx'
myx, myy, slope, intercept = readxlsx(filename)
popt, pcov = curve_fit(linear, myx, myy, p0 = (1.0, 1.0))
myslopenoisy = popt[0]
myinterceptnoisy = popt[1]
myfittedynoisy = []
for i in range(len(myx)):
	myfittedynoisy.append(myslopenoisy*myx[i]+myinterceptnoisy)

plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.plot(myx, myy,linestyle=':')
plt.plot(myx, myfittedynoisy)
plt.savefig('noisyfigure.png')
plt.clf()




@

Now execute the dolinear script 

\begin{pycode}
from dolinear import *
\end{pycode}

You can access the results/variables from the dolinear script

\begin{figure}
\includegraphics[scale=0.75]{figure.png}
\caption{No noise added}
\end{figure}

\begin{figure}
\includegraphics[scale=0.75]{noisyfigure.png}
\caption{Noise added}
\end{figure}

The equation was $y = slope*x + intercept$.  The slope value used
was \py{slope} and intercept value was \py{intercept}.

if noise was NOT added, the 
calculated slope was \py{myslope} and the 
calculated intercept was \py{myintercept}

When noise was added to the $y$ and when curve fit was used, the 
calculated slope was \py{myslopenoisy} and the 
calculated intercept was \py{myinterceptnoisy}


\end{document}

