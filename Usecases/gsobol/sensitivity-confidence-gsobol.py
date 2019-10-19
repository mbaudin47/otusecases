# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 11:48:42 2017

@author: c61372

Analyse la précision des estimateurs des intervalles de confiance 
pour la fonction G-Sobol
"""

#! /usr/bin/env python

from __future__ import print_function
import openturns as ot
from gsobollib import (
        gsobolSAExact, 
        gsobolDistribution, gsobol
)
import numpy as np
import pylab as pl
import openturns.viewer as otv

a = np.array([0,9,99])
d = len(a)

# Distribution uniforme associée au cas-test GSobol
distribution = gsobolDistribution(d)

# Indices de sensibilité exacts
[muexact,vexact,sexact,stexact] = gsobolSAExact(a)

# Taille du plan d'expérience de base pour estimer S et ST
sampleSize = 1000

# Nombre de répétition de l'expérience
nrepetitions = 500

# Estimations des indices du premier ordre
sampleFirstMartinez = ot.Sample(nrepetitions,3)

# Estimations des indices totaux
sampleTotalMartinez = ot.Sample(nrepetitions,3)
for i in range(nrepetitions):
    inputDesign = ot.SobolIndicesExperiment(distribution, sampleSize).generate()
    outputDesign = gsobol(inputDesign,a)
    sensitivity_algorithm = ot.SaltelliSensitivityAlgorithm(
        inputDesign, outputDesign, sampleSize)
    fo = sensitivity_algorithm.getFirstOrderIndices()
    to = sensitivity_algorithm.getTotalOrderIndices()
    for j in range(d):
        sampleFirstMartinez[i,j] = fo[j]
    for j in range(d):
        sampleTotalMartinez[i,j] = to[j]

fig = pl.figure(figsize=(12, 8))
for j in range(d):
    ax = fig.add_subplot(2, 3, 1+j)
    graph = ot.HistogramFactory().build(sampleFirstMartinez[:,j]).drawPDF()
    graph.setXTitle("S%d" % (d))
    graph.setLegends([""])
    _ = otv.View(graph, figure=fig, axes=[ax])
    ax = fig.add_subplot(2,3,4+j)
    graph = ot.HistogramFactory().build(sampleTotalMartinez[:,j]).drawPDF()
    graph.setXTitle("ST%d" % (d))
    graph.setLegends([""])
    _ = otv.View(graph, figure=fig, axes=[ax])
_ = fig.suptitle("Martinez - N=%d - Repetitions = %d" % (sampleSize,nrepetitions))

# Récupère l'intervalle de confiance bootstrap pour le dernier échantillon
alpha = sensitivity_algorithm.getConfidenceLevel()
foInterval = sensitivity_algorithm.getFirstOrderIndicesInterval()
foIntervalMin = foInterval.getLowerBound()
foIntervalMax = foInterval.getUpperBound()
toInterval = sensitivity_algorithm.getTotalOrderIndicesInterval()
toIntervalMin = toInterval.getLowerBound()
toIntervalMax = toInterval.getUpperBound()

# Compare les intervalles bootstrap pour le dernier échantillon 
# et les quantiles issus des répétitions
for j in range(d):
    # Calcule les quantiles empiriques
    sampleFirst = sampleFirstMartinez[:,j]
    foMinj = sampleFirst.computeQuantile((1-alpha)/2)[0]
    foMaxj = sampleFirst.computeQuantile(1-(1-alpha)/2)[0]
    sampleTotal = sampleTotalMartinez[:,j]
    toMinj = sampleTotal.computeQuantile((1-alpha)/2)[0]
    toMaxj = sampleTotal.computeQuantile(1-(1-alpha)/2)[0]
    print("X%d" % (j))
    print("   First, Bootstrap=[%.4f,%.4f], Sample=[%.4f,%.4f]" % (foIntervalMin[j],foIntervalMax[j],foMinj,foMaxj))
    print("   Total, Bootstrap=[%.4f,%.4f], Sample=[%.4f,%.4f]" % (toIntervalMin[j],toIntervalMax[j],toMinj,toMaxj))

fig = pl.figure(figsize=(12, 8))
for j in range(d):
    # First order
    ax = fig.add_subplot(2, 3, 1+j)
    sampleJ = sampleFirstMartinez[:,j]
    graph = ot.HistogramFactory().build(sampleJ).drawPDF()
    graph.setXTitle("S%d" % (d))
    mu = sexact[j] # Valeur exacte
    # TODO : mettre la valeur issue de l'estimateur asymptotique
    sigma = sampleJ.computeStandardDeviationPerComponent()[0] 
    distribution = ot.Normal(mu,sigma)
    graphPDF = distribution.drawPDF()
    graphPDF.setColors(["blue"])
    graph.add(graphPDF)
    graph.setLegends([""])
    _ = otv.View(graph, figure=fig, axes=[ax])
    # Total order
    sampleJ = sampleTotalMartinez[:,j]
    ax = fig.add_subplot(2,3,4+j)
    graph = ot.HistogramFactory().build(sampleJ).drawPDF()
    mu = stexact[j] # Valeur exacte
    # TODO : mettre la valeur issue de l'estimateur asymptotique
    sigma = sampleJ.computeStandardDeviationPerComponent()[0] 
    distribution = ot.Normal(mu,sigma)
    graphPDF = distribution.drawPDF()
    graphPDF.setColors(["blue"])
    graph.setXTitle("ST%d" % (d))
    graph.add(graphPDF)
    graph.setLegends([""])
    _ = otv.View(graph, figure=fig, axes=[ax])
_ = fig.suptitle("Martinez - N=%d - Repetitions = %d" % (sampleSize,nrepetitions))
