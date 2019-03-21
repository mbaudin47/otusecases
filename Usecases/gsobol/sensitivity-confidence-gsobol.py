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
from numpy import zeros, sqrt, array, linspace
from pylab import plot, show, xlabel, ylabel, xscale, yscale, legend, title, savefig, subplot, hist

a = array([0,9,99])
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
sampleFirstMartinez = zeros((nrepetitions,3))

# Estimations des indices totaux
sampleTotalMartinez = zeros((nrepetitions,3))
for i in range(nrepetitions):
    inputDesign = ot.SaltelliSensitivityAlgorithm.Generate(
        distribution, int(sampleSize))
    outputDesign = gsobol(inputDesign,a)
    sensitivity_algorithm = ot.SaltelliSensitivityAlgorithm(
        inputDesign, outputDesign, int(sampleSize))
    fo = sensitivity_algorithm.getFirstOrderIndices()
    to = sensitivity_algorithm.getTotalOrderIndices()
    for j in range(d):
        sampleFirstMartinez[i,j] = fo[j]
    for j in range(d):
        sampleTotalMartinez[i,j] = to[j]

subplot(2,3,1)
title("Martinez - N=%d - Repetitions = %d" % (sampleSize,nrepetitions))
for j in range(d):
    subplot(2,3,1+j)
    hist(sampleFirstMartinez[:,j],histtype="step",normed=True)
    xlabel("S%d" % (d))
    ylabel("Density")
    subplot(2,3,4+j)
    hist(sampleTotalMartinez[:,j],histtype="step",normed=True)
    xlabel("ST%d" % (d))
    ylabel("Density")
show()

# Récupère l'intervalle de confiance bootstrap pour le dernier échantillon
alpha = sensitivity_algorithm.getBootstrapConfidenceLevel()
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
    sampleFirst = ot.Sample(sampleFirstMartinez[:,j],1)
    foMinj = sampleFirst.computeQuantile((1-alpha)/2)[0]
    foMaxj = sampleFirst.computeQuantile(1-(1-alpha)/2)[0]
    sampleTotal = ot.Sample(sampleFirstMartinez[:,j],1)
    toMinj = sampleTotal.computeQuantile((1-alpha)/2)[0]
    toMaxj = sampleTotal.computeQuantile(1-(1-alpha)/2)[0]
    print("X%d" % (j))
    print("   First, Bootstrap=[%.4f,%.4f], Sample=[%.4f,%.4f]" % (foIntervalMin[j],foIntervalMax[j],foMinj,foMaxj))
    print("   Total, Bootstrap=[%.4f,%.4f], Sample=[%.4f,%.4f]" % (toIntervalMin[j],toIntervalMax[j],toMinj,toMaxj))

# Pour chaque estimateur, compare la répartition empirique et la loi exacte
subplot(2,3,1)
title("Martinez - N=%d - Repetitions = %d" % (sampleSize,nrepetitions))
for j in range(d):
    # Indice du premier ordre
    subplot(2,3,1+j)
    sampleJ = sampleFirstMartinez[:,j]
    hist(sampleJ,histtype="step",normed=True)
    sampleFirst = ot.Sample(sampleJ,1)
    mu = sexact[j] # Valeur exacte
    # TODO : mettre la valeur issue de l'estimateur asymptotique
    sigma = sampleFirst.computeStandardDeviationPerComponent()[0] 
    loiFo = ot.Normal(mu,sigma)
    x = linspace(mu-3*sigma,mu+3*sigma,50)
    y = loiFo.computePDF(ot.Sample(x,1))
    plot(x,y,"r-")
    xlabel("S%d" % (j))
    ylabel("Density")
    # Indice du total
    subplot(2,3,4+j)
    sampleJ = sampleTotalMartinez[:,j]
    hist(sampleJ,histtype="step",normed=True)
    sampleFirst = ot.Sample(sampleJ,1)
    mu = stexact[j] # Valeur exacte
    # TODO : mettre la valeur issue de l'estimateur asymptotique
    sigma = sampleFirst.computeStandardDeviationPerComponent()[0] 
    loiFo = ot.Normal(mu,sigma)
    x = linspace(mu-3*sigma,mu+3*sigma,50)
    y = loiFo.computePDF(ot.Sample(x,1))
    plot(x,y,"r-")
    xlabel("ST%d" % (j))
    ylabel("Density")
show()
