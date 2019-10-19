#!/usr/bin/env python
# coding: utf-8

import openturns as ot
from openturns.viewer import View

# 1. Define the random variables
myParam = ot.GumbelAB(1013., 558.)
Q = ot.ParametrizedDistribution(myParam)
otLOW = ot.TruncatedDistribution.LOWER
Q = ot.TruncatedDistribution(Q, 0, otLOW)
Ks = ot.Normal(30.0, 7.5)
Ks = ot.TruncatedDistribution(Ks, 0, otLOW)

# 2. Define the function
formulas = ['(Q/(Ks*300.*sqrt(0.001)))^(3./5.)']
g = ot.SymbolicFunction(['Q','Ks'],formulas)

# 3. Create the joint distribution
inputDistribution = ot.ComposedDistribution((Q,Ks))
inputRandomVector = ot.RandomVector(inputDistribution)
outputRandomVector = ot.CompositeRandomVector(g, inputRandomVector)

# 4. Simple Monte-Carlo sampling
samplesize=500
outputSample=outputRandomVector.getSample(samplesize)

# 6. Plot the histogram
histoGraph = ot.HistogramFactory().build(outputSample).drawPDF()
histoGraph.setTitle("Histogramme de la hauteur")
histoGraph.setXTitle("H (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setLegends([""])
View(histoGraph)
