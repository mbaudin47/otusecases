#!/usr/bin/env python
# coding: utf-8

import openturns as ot
from math import sqrt
from openturns.viewer import View

# 1. The function G
def functionCrue8vars(X) :
    Q, Ks, Zv, Zm, Hd, Zb, L, B = X
    Zd = Zb + Hd
    alpha = (Zm - Zv)/L
    H = (Q/(Ks*B*sqrt(alpha)))**(3.0/5.0)
    Zc = H + Zv
    S = Zc - Zd
    return [S]

def functionCrue4VarsStochastic(X):
    Q, Ks, Zv, Zm = X
    # 1. Creation of the problem function
    f8v = ot.PythonFunction(8, 1, functionCrue8vars)
    g = ot.ParametricFunction(f8v, [0, 1, 2, 3], X)
    # 2. Random vector definition
    Hd = ot.Uniform(4.,14.)
    Zb = ot.Uniform(50.,60.)
    L = ot.Uniform(1000.,10000.)
    B = ot.Uniform(50.,500.)
    inputvector = ot.ComposedDistribution([Hd, Zb, L, B])
    inputRV = ot.RandomVector(inputvector)
    S = ot.CompositeRandomVector(g, inputRV)
    # 3. Sample output
    sampleSize = 10
    outputSample = S.getSample(sampleSize)
    Smean = outputSample.computeMean()[0]
    return [Smean]

# Creation of the problem function
g = ot.PythonFunction(4, 1, functionCrue4VarsStochastic) 

# See the stochastic code in action
# Two consecutive calls do not produce the same result
X = [1013.,30.,50.,55.]
print(g(X))
print(g(X))
'''
Do *not* configure a cache !

g = ot.MemoizeFunction(g)

Otherwise, two calls will alwas produce the same result.
'''

# 2. Random vector definition
myParam = ot.GumbelAB(1013., 558.)
Q = ot.ParametrizedDistribution(myParam)
otLOW = ot.TruncatedDistribution.LOWER
Q = ot.TruncatedDistribution(Q, 0, otLOW)
Ks = ot.Normal(30.0, 7.5)
Ks = ot.TruncatedDistribution(Ks, 0, otLOW)
Zv = ot.Uniform(49.0, 51.0)
Zm = ot.Uniform(54.0, 56.0)

# 3. Create the joint distribution
inputDistribution = ot.ComposedDistribution((Q,Ks,Zv,Zm))
inputRandomVector = ot.RandomVector(inputDistribution)
outputRandomVector = ot.CompositeRandomVector(g, inputRandomVector)

# 4. Simple Monte-Carlo sampling
samplesize=500
outputSample=outputRandomVector.getSample(samplesize)

# 6. Plot the histogram
histoGraph = ot.HistogramFactory().build(outputSample).drawPDF()
histoGraph.setTitle("Histogramme de la surverse")
histoGraph.setXTitle("S (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setLegends([""])
View(histoGraph)

