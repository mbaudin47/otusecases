# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:28:26 2016

Objectif :
Montrer comment estimer les indices avec un plan Monte-Carlo simple.
Voir la différence avec un plan de type séquence à faible discrépance.
"""

import openturns as ot
from math import pi
dim = 3
size = 20
fla = "sin(X1) + 7*sin(X2)^2 + 0.1*X3^4*sin (X1)"
g = ot.SymbolicFunction ([ "X1", "X2", "X3"], [fla ])
X = ot.ComposedDistribution ([ ot.Uniform (-pi , pi)] * dim )

alpha = 0.05 # i.e. 90% confidence interval
epsilon = 0.2 # Confidence interval length
blocksize = 50 # size of Sobol experiment at each iteration
batchsize = 16 # number of points evaluated simultaneously

estimator = ot.SaltelliSensitivityAlgorithm()
estimator.setUseAsymptoticDistribution(True)
algo = ot.SobolSimulationAlgorithm(X, g, estimator)
algo.setMaximumOuterSampling(100) # number of iterations
algo.setBlockSize(blocksize) 
algo.setBatchSize(batchsize) 
algo.setIndexQuantileLevel(alpha) # alpha
algo.setIndexQuantileEpsilon(epsilon) # epsilon
algo.run()

result = algo.getResult()
fo = result.getFirstOrderIndicesEstimate()
to = result.getTotalOrderIndicesEstimate()
print("First order = %s" % (fo))
print("Total order = %s" % (to))
