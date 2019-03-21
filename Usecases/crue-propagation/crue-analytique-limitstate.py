# -*- coding: utf-8 -*-
"""
Scénario : définir le code G - Analytique
Fonction G : Analytique
PDF des 4 variables d'entrée X
Matrice de scatter plot de l'échantillon de l'entrée X
Estimation d'une probabilité de dépassement de seuil, 
  et son intervalle de confiance à 95%.
Algorithme : MonteCarlo
Histogramme de la surverse
"""

from openturns.viewer import View
from numpy import inf
from openturns import (
    ComposedDistribution, Gumbel, 
    RandomVector, Normal, NumericalMathFunction, MonteCarlo, 
    TruncatedDistribution, Event, Uniform, GreaterOrEqual
)
import time
from openturns import VisualTest
from math import ceil, log
from pylab import plot, show, xscale, xlabel, ylabel, legend, title, grid
from numpy import array, sqrt, where

# 1. The function G

inputs=["Q","Ks","Zv","Zm"]
outputs=["S"]
formula=["Zv+(Q/(300*Ks*sqrt((Zm-Zv)/5000)))^0.6-3-55.5"]
f = NumericalMathFunction(inputs,outputs,formula)
f.enableHistory()

# 2. Random vector definition
Q = Gumbel(1./558., 1013.)
Q = TruncatedDistribution(Q, 0, inf)
K_s = Normal(30.0, 7.5)
K_s = TruncatedDistribution(K_s, 0, inf)
Z_v = Uniform(49.0, 51.0)
Z_m = Uniform(54.0, 56.0)

# 4. Create the joint distribution function, 
#    the output and the event. 
inputRandomVector = ComposedDistribution([Q, K_s, Z_v, Z_m])
outputRandomVector = RandomVector(f, RandomVector(inputRandomVector))

# Define event
threshold = -2.
eventF = Event(outputRandomVector, GreaterOrEqual(), threshold) 


# 5. Create the Monte-Carlo algorithm
algoProb = MonteCarlo(eventF)
montecarlosize=1000000
algoProb.setMaximumOuterSampling(montecarlosize)
algoProb.setMaximumCoefficientOfVariation(0.1)

callstart_time=time.time()
algoProb.run()
callstop_time=time.time()
elapsed_time=callstop_time - callstart_time
print "Elapsed time: %.2f (s)" % (elapsed_time)

# 6. Get the results
resultAlgo = algoProb.getResult()
Neval = f.getEvaluationCallsNumber()
print "Number of function calls =", Neval
Pf = resultAlgo.getProbabilityEstimate()
print "Failure Probability = %e" % (Pf)
c95 = resultAlgo.getConfidenceLength()
pmin=Pf-0.5*c95
pmax=Pf+0.5*c95
print "within [%e,%e] at 95%%" % (pmin,pmax)
cov=resultAlgo.getCoefficientOfVariation()
print "Coeff. of Var. = %.2f" % (cov)



# 8. Plot the histogram - better
def computeHistogramBoundsForLimitState(threshold, outComputedPoints):
    # Returns the bounds for the histogram printed 
    # for the Limit State analysis
    xmin=outComputedPoints.computeQuantile(0.01)[0]
    xmax=outComputedPoints.computeQuantile(0.99)[0]
    delta=xmax-xmin
    if threshold> xmax:
        xmax=threshold+0.2*delta
    if threshold< xmin:
        xmin=threshold+0.2*delta
    return xmin,xmax

def histogramClassNumberSturges(n):
    # Returns Sturges rule for number of classes 
    # of an histogram
    nclasses=int(ceil(log(n,2)+1))
    return nclasses

def histogramClassNumberScott(n,sigma,xmin,xmax):
    # Returns Scott's rule for number of classes
    # of an histogram
    h=3.5*sigma/n**(1./3.)
    nclasses=int(ceil((xmax-xmin)/h))
    return nclasses


# 8.1 Prepare data
outComputedPoints = f.getOutputHistory().getSample()
n=outComputedPoints.getSize()
sigma=outComputedPoints.computeStandardDeviationPerComponent()[0]
xmin=outComputedPoints.getMin()[0]
xmax=outComputedPoints.getMax()[0]
barsNumber = histogramClassNumberScott(n,sigma,xmin,xmax)

# 8.3 Plot the histogram
histoGraph = VisualTest.DrawHistogram(outComputedPoints,barsNumber)
histoGraph.setTitle("Histogramme de la surverse")
histoGraph.setXTitle("S (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setLegends([""])
boundingbox=histoGraph.getBoundingBox()
xmin,xmax=computeHistogramBoundsForLimitState(threshold, outComputedPoints)
boundingbox[0]=xmin
boundingbox[1]=xmax
histoGraph.setBoundingBox(boundingbox)
View(histoGraph).show()

# 8.2 Plot a red line at threshold
drawable0=histoGraph.getDrawables()[0]
bbox=drawable0.getBoundingBox()
ymax=bbox[3]
plot([threshold,threshold],[0,0.1*ymax],"r-",linewidth=2)
show()

# 9. Graphique MonteCarlo - échelle log

# 9.1 Get data
algonconvergence=algoProb.getConvergenceStrategy()
echantillonEstimateur = algonconvergence.getSample()
N=echantillonEstimateur.getSize()
data=array(echantillonEstimateur[:,0])
datavar=array(echantillonEstimateur[:,1])

# Extract the data which is consistent with log-scale
datamin=data-sqrt(datavar) # Lower bound
datamax=data+sqrt(datavar) # Upper bound
first=where(data-sqrt(datavar)>0)[0][0] # First positive lower bound
iterationsindices=range(first,N)

# 9.2 Graphics
xscale("log")
xlabel("Number of iterations")
ylabel("Estimate")
plot(iterationsindices,data[first:],"r-", label="Probability Estimate")
plot(iterationsindices,datamax[first:],"g-", label="Bounds")
plot(iterationsindices,datamin[first:],"g-")
title("MonteCarlo convergence graph at level 0.95")
legend()
grid(True)
show()
