# -*- coding: utf-8 -*-
"""
Created on Thu Jun 01 11:48:42 2017

@author: c61372

Analyse la convergence des estimateurs des indices de sensibilité 
pour la fonction G-Sobol.
"""

#! /usr/bin/env python

from __future__ import print_function
import openturns as ot
from gsobollib import (
        gsobolSAExact, 
        gsobolDistribution, gsobol
)
from numpy import zeros, sqrt, array
from pylab import plot, show, xlabel, ylabel, xscale, yscale, legend, title, savefig

def computeAbsoluteError(sensitivity_algorithm,a):
    # Get approximate indices
    fo = sensitivity_algorithm.getFirstOrderIndices()
    to = sensitivity_algorithm.getTotalOrderIndices()
    # Get exact indices
    [muexact,vexact,sexact,stexact] = gsobolSAExact(a)
    # Compute absolute errors
    AES1 = abs(fo[0]-sexact[0])
    AES2 = abs(fo[1]-sexact[1])
    AES3 = abs(fo[2]-sexact[2])
    AEST1 = abs(to[0]-stexact[0])
    AEST2 = abs(to[1]-stexact[1])
    AEST3 = abs(to[2]-stexact[2])
    absErrFirst = max([AES1,AES2,AES3])
    absErrTotal = max([AEST1,AEST2,AEST3])
    return [absErrFirst,absErrTotal]

def plotAbsoluteError(algorithmName,sampleSize,absErrorFirst,absErrorTotal):
    title(algorithmName)
    plot(sampleSize,1./sqrt(sampleSize),"-", label="1/sqrt(n)")
    plot(sampleSize,absErrorFirst,"o", label="First order")
    plot(sampleSize,absErrorTotal,"o", label="Total order")
    xlabel("N")
    ylabel("Absolute error")
    xscale("log")
    yscale("log")
    legend()
    show()
    return None

a = array([0,9,99])
nx = len(a)

distribution = gsobolDistribution(nx)

# Size of simulation
nloops = 15
sampleSize = zeros((nloops,1))
sampleSize[0] = 10
for i in range(nloops):
    if (i>0):
        sampleSize[i] = 2*sampleSize[i-1]

#########################################
# 1. Saltelli

absErrorFirst = zeros((nloops,1))
absErrorTotal = zeros((nloops,1))
for i in range(nloops):
    inputDesign = ot.SaltelliSensitivityAlgorithm.Generate(
        distribution, int(sampleSize[i]))
    outputDesign = gsobol(inputDesign,a)
    sensitivity_algorithm = ot.SaltelliSensitivityAlgorithm(
        inputDesign, outputDesign, int(sampleSize[i]))
    # Compute accuracy
    absErrFirst_i,absErrTotal_i = computeAbsoluteError(sensitivity_algorithm,a)
    absErrorFirst[i] = absErrFirst_i
    absErrorTotal[i] = absErrTotal_i

plotAbsoluteError("Gsobol-Saltelli",sampleSize,absErrorFirst,absErrorTotal)

#########################################
# 2. Martinez

absErrorFirst = zeros((nloops,1))
absErrorTotal = zeros((nloops,1))
for i in range(nloops):
    computeSecondOrder = False
    inputDesign = ot.MartinezSensitivityAlgorithm.Generate(
        distribution, int(sampleSize[i]),computeSecondOrder)
    outputDesign = gsobol(inputDesign,a)
    sensitivity_algorithm = ot.MartinezSensitivityAlgorithm(
        inputDesign, outputDesign, int(sampleSize[i]))
    # Compute accuracy
    absErrFirst_i,absErrTotal_i = computeAbsoluteError(sensitivity_algorithm,a)
    absErrorFirst[i] = absErrFirst_i
    absErrorTotal[i] = absErrTotal_i

plotAbsoluteError("Gsobol-Martinez",sampleSize,absErrorFirst,absErrorTotal)

#########################################
# 3. Jansen

absErrorFirst = zeros((nloops,1))
absErrorTotal = zeros((nloops,1))
for i in range(nloops):
    computeSecondOrder = False
    inputDesign = ot.JansenSensitivityAlgorithm.Generate(
        distribution, int(sampleSize[i]),computeSecondOrder)
    outputDesign = gsobol(inputDesign,a)
    sensitivity_algorithm = ot.JansenSensitivityAlgorithm(
        inputDesign, outputDesign, int(sampleSize[i]))
    # Compute accuracy
    absErrFirst_i,absErrTotal_i = computeAbsoluteError(sensitivity_algorithm,a)
    absErrorFirst[i] = absErrFirst_i
    absErrorTotal[i] = absErrTotal_i

plotAbsoluteError("Gsobol-Jansen",sampleSize,absErrorFirst,absErrorTotal)

#########################################
# 4. MauntzKucherenko

absErrorFirst = zeros((nloops,1))
absErrorTotal = zeros((nloops,1))
for i in range(nloops):
    computeSecondOrder = False
    inputDesign = ot.MauntzKucherenkoSensitivityAlgorithm.Generate(
        distribution, int(sampleSize[i]),computeSecondOrder)
    outputDesign = gsobol(inputDesign,a)
    sensitivity_algorithm = ot.MauntzKucherenkoSensitivityAlgorithm(
        inputDesign, outputDesign, int(sampleSize[i]))
    # Compute accuracy
    absErrFirst_i,absErrTotal_i = computeAbsoluteError(sensitivity_algorithm,a)
    absErrorFirst[i] = absErrFirst_i
    absErrorTotal[i] = absErrTotal_i

plotAbsoluteError("Gsobol-MauntzKucherenko",sampleSize,absErrorFirst,absErrorTotal)
