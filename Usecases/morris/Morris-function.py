#!/usr/bin/env python
# coding: utf-8

# # The Morris function
# 
# ## Reference
# 
# * M. D. Morris, 1991, Factorial sampling plans for preliminary computational experiments,Technometrics, 33, 161--174.
# 

import openturns as ot
import otmorris
from openturns.viewer import View


# Setting the seed allows to produce always the same Morris parameters.

ot.RandomGenerator.SetSeed(1)

alpha = ot.DistFunc.rNormal(10)
beta = ot.DistFunc.rNormal(14*6)
gamma = ot.DistFunc.rNormal(20*14)
b0 = ot.DistFunc.rNormal()

g = ot.Function(otmorris.MorrisFunction(alpha, beta, gamma, b0))

U = [ot.Uniform(0,1)]*20

X = ot.ComposedDistribution(U)

n = 500

sampleX = X.getSample(n)

sampleY = g(sampleX)

graph = ot.VisualTest_DrawHistogram(sampleY)
graph.setLegends([""])
graph.setTitle("Morris function")
graph.setXTitle("g(X)")
graph.setYTitle("PDF")
View(graph)

