#!/usr/bin/env python
# coding: utf-8

import openturns as ot

R = ot.Normal(4., 1.)
R.setDescription("R")

S = ot.Normal(2., 1.)
S.setDescription("S")

g = ot.SymbolicFunction(["R","S"],["R-S"])

inputvector = ot.ComposedDistribution([R,S])
inputRV = ot.RandomVector(inputvector)
outputRV = ot.CompositeRandomVector(g, inputRV)
eventF = ot.Event(outputRV, ot.GreaterOrEqual(), 0) 

# Create the Monte-Carlo algorithm
algoProb = ot.ProbabilitySimulationAlgorithm(eventF)
algoProb.setMaximumOuterSampling(1000)
algoProb.setMaximumCoefficientOfVariation(0.01)
algoProb.run()

# Get the results
resultAlgo = algoProb.getResult()
neval = g.getEvaluationCallsNumber()
print("Number of function calls = %d" %(neval))
pf = resultAlgo.getProbabilityEstimate()
print("Failure Probability = %.4f" % (pf))
level = 0.95
c95 = resultAlgo.getConfidenceLength(level)
pmin=pf-0.5*c95
pmax=pf+0.5*c95
print("%.1f %% confidence interval :[%.4f,%.4f] " % (level*100,pmin,pmax))
