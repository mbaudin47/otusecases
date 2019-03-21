# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:28:26 2016

Objectif :
Faire une étude en tendance centrale sur la formule de 
la fleche d'un tube.
"""

from openturns import (
NumericalMathFunction, ComposedDistribution, 
Uniform, Normal, Triangular, RandomVector, 
SensitivityAnalysis, CorrelationMatrix, NormalCopula, 
Graph, Pairs, IndependentCopula
)
from centraldispersion import (
centralDispersionByMonteCarlo, 
centralDispersionPrintResults, 
centralDispersionPrintParameters
)
from openturns.viewer import View
from openturns import VisualTest
from openturns import CorrelationAnalysis_SRC
from plotsa import plotSRCindices, plotSobolIndices

# 1. Define the function
dim = 6
inputsvars=["F","L","a","D","d","E"]
formula = ["-F*a^2*(L-a)^2/(3*E*L*3.14159*(D^4-d^4)/32)"]
model = NumericalMathFunction (inputsvars,formula)

# 2. Define the probabilistic model
XF=Normal(1,0.1)
XL=Normal(1.5,0.01)
Xa=Uniform(0.7,1.2)
XD=Triangular(0.75,0.8,0.85)
Xd=Triangular(0.09,0.1,0.11)
XE=Normal(200000,2000)
XF.setDescription(["Sollicitation"])
XL.setDescription(["Longueur"])
Xa.setDescription(["Appui"])
XD.setDescription(["Diametre externe"])
Xd.setDescription(["Diametre interne"])
XE.setDescription(["Module d'Young"])

if False:
    View(XF.drawPDF()).show()
    View(XL.drawPDF()).show()
    View(Xa.drawPDF()).show()
    View(XE.drawPDF()).show()
    View(XD.drawPDF()).show()
    View(Xd.drawPDF()).show()

# Create the Spearman correlation matrix of the input random vector
if False:
    RS = CorrelationMatrix(6)
    RS[2,5] = 0.8
    # Evaluate the correlation matrix of the Normal copula from RS
    R = NormalCopula.GetCorrelationFromSpearmanCorrelation(RS)
    # Create the Normal copula parametrized by R
    mycopula = NormalCopula(R)
else:
    mycopula=IndependentCopula(6)

# 3. Create the joint distribution function, 
#    the output and the event. 
inputDistribution = ComposedDistribution([XF,XL,Xa,XD,Xd,XE], mycopula)

# Create the input random vector
inputRandomVector = RandomVector(inputDistribution)

outputRandomVector = RandomVector(model, inputRandomVector)
outputDim = model.getOutputDimension()

# 3.1 Draw scatter plot matrix
if False:
    samplesize=100
    mytitle='Sample - N=%d' % (samplesize)
    sample = inputRandomVector.getSample(samplesize)
    myGraph = Graph(mytitle, ' ', ' ', True, '')
    # Create the Pairs
    mydescription=sample.getDescription()
    myPairs = Pairs(sample, mytitle, mydescription, 'blue', 'bullet')
    myGraph.add(myPairs)
    View(myGraph).show()

# 4. Estimate central dispersion

# User parameters
blocksize=1000
maxcov=0.001 # Criteria B
maxcalls=100000 # Criteria A
maxelapsetime=5 # Criteria C
alpha=0.05 # Niveau de confiance de l'intervalle

centralDispersionPrintParameters(blocksize, maxcov, maxcalls,maxelapsetime,alpha)

print ""
print "Use centraldispersion lib"

outputSample, criteria=centralDispersionByMonteCarlo(outputRandomVector, blocksize,maxcov,maxcalls,maxelapsetime)
if (criteria==1):
    print "Reached number of calls"
elif (criteria==2):
    print "Reached required precision"
elif (criteria==3):
    print "Reached maximum elapsed time"
    
centralDispersionPrintResults(outputSample,alpha)

# 5. Analyse central dispersion

# Get the empirical mean and standard deviations
empiricalMean = outputSample.computeMean()
empiricalSd = outputSample.computeStandardDeviationPerComponent()
empiricalMin=outputSample.getMin()
empiricalMax=outputSample.getMax()
empiricalSkew=outputSample.computeSkewness()
for i in range(outputDim):
    print "Sortie #%d" % (i)
    print "\t Mean=%e" % (empiricalMean[i])
    print "\t Sd.Dev.=%e" % (empiricalSd[i])
    print "\t Min=%e" % (empiricalMin[i])
    print "\t Max=%e" % (empiricalMax[i])
    print "\t Skewness=%e" % (empiricalSkew[i])

# 7. Plot the histogram
if False:
    histoGraph = VisualTest.DrawHistogram(outputSample)
    samplesize=outputSample.getSize()
    histoGraph.setTitle("Histogramme de la flèche - N=%d" % (samplesize))
    histoGraph.setXTitle("Flèche")
    histoGraph.setYTitle("Fréquence")
    histoGraph.setLegends([""])
    View(histoGraph).show()

# CDF
if False:
    empiricalDFGraph = VisualTest.DrawEmpiricalCDF(outputSample,empiricalMin[0],empiricalMax[0])
    samplesize=outputSample.getSize()
    empiricalDFGraph.setTitle("Histogramme de la flèche - N=%d" % (samplesize))
    empiricalDFGraph.setXTitle("Flèche")
    empiricalDFGraph.setYTitle("Fréquence Cumulée")
    empiricalDFGraph.setLegends([""])
    View(empiricalDFGraph).show()

# 8. Compute the Standard Regression Coefficients 
# with respect to the inputs.
print "Standard Regression Coefficients:"
montecarlosize = 1000
print "Nombre de simulations : %d" % (montecarlosize)
inputSample=inputRandomVector.getSample(montecarlosize)
outputSample = model(inputSample)
srcindices = CorrelationAnalysis_SRC(inputSample,outputSample)
ylabdescr=outputRandomVector.getDescription()[0]
print "Y=", ylabdescr
inputDim=model.getInputDimension()
print "SRC indices:"
for i in range(inputDim):
    descr=inputDistribution.getDescription()[i]
    print "X#%d, %s, SRC=%.2f %%" % (i,descr,srcindices[i]*100)

print "Sum:%.2f %%" % (sum(srcindices)*100)

# 9. Plot SRC indices
plotSRCindices(srcindices,inputDistribution)


# 10. Compute Sobol indices
montecarlosize = 10000
inputSample1 = inputDistribution.getSample(montecarlosize)
inputSample2 = inputDistribution.getSample(montecarlosize)
print "Nombre de simulations : %d" % (montecarlosize)
sensitivityAnalysis = SensitivityAnalysis(inputSample1, inputSample2, model)
firstOrderIndices = sensitivityAnalysis.getFirstOrderIndices()
totalOrderIndices = sensitivityAnalysis.getTotalOrderIndices()
ylabdescr=outputRandomVector.getDescription()[0]
print "Y=", ylabdescr

# 8. Print some indices
inputDim=model.getInputDimension()
print "First order Sobol indices:"
for i in range(inputDim):
    descr=inputDistribution.getDescription()[i]
    print "X#%d, %s, S=%.2f %%" % (i,descr,firstOrderIndices[i]*100)

print "Sum:%.2f %%" % (sum(firstOrderIndices)*100)

print "Total order Sobol indices:"
for i in range(inputDim):
    descr=inputDistribution.getDescription()[i]
    print "X#%d, %s, ST=%.2f %%" % (i,descr,totalOrderIndices[i]*100)

print "Sum:%.2f %%" % (sum(totalOrderIndices)*100)

# 9. Plot SRC indices
plotSobolIndices(firstOrderIndices,totalOrderIndices,inputDistribution)
