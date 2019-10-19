#!/usr/bin/env python
# coding: utf-8

import openturns as ot
from openturns.viewer import View

# 2. Create the Input and Output random variables
myParam = ot.GumbelAB(1013., 558.)
QGumbel = ot.ParametrizedDistribution(myParam)
Q = ot.TruncatedDistribution(QGumbel, 0, ot.TruncatedDistribution.LOWER)
KsNormal = ot.Normal(30.0, 7.5)
Ks = ot.TruncatedDistribution(KsNormal, 0, ot.TruncatedDistribution.LOWER)
Zv = ot.Uniform(49.0, 51.0)
Zm = ot.Uniform(54.0, 56.0)
#
Hd = ot.Uniform(7., 9.) # Hd = 3.0;
Zb = ot.Triangular(55.0, 55.5, 56.0) # Zb = 55.5
L = ot.Triangular(4990, 5000., 5010.) # L = 5.0e3;
B = ot.Triangular(295., 300., 305.) # B = 300.0

Q.setDescription(["Q (m3/s)"])
Ks.setDescription(["Ks (m^(1/3)/s)"])
Zv.setDescription(["Zv (m)"])
Zm.setDescription(["Zm (m)"])
Hd.setDescription(["Hd (m)"])
Zb.setDescription(["Zb (m)"])
L.setDescription(["L (m)"])
B.setDescription(["B (m)"])


formulas = ['(Q/(Ks*B*sqrt((Zm-Zv)/L)))^(3.0/5.0)+Zv-Zb-Hd']
g = ot.SymbolicFunction(['Q','Ks','Zv','Zm','Hd','Zb','L','B'],formulas)

# 3. Create the joint distribution
inputDistribution = ot.ComposedDistribution((Q,Ks,Zv,Zm,Hd,Zb,L,B))
inputRandomVector = ot.RandomVector(inputDistribution)
outputRandomVector = ot.CompositeRandomVector(g, inputRandomVector)

# 4. Simple Monte-Carlo sampling
samplesize=500
outputSample=outputRandomVector.getSample(samplesize)

# 6. Plot the histogram
histoGraph = ot.HistogramFactory().build(outputSample[:,0]).drawPDF()
histoGraph.setTitle("Histogramme de la surverse")
histoGraph.setXTitle("S (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setLegends([""])
View(histoGraph)

