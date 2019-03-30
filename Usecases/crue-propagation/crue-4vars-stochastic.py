#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import openturns as ot
import otguibase

Study_Crue4Stochastic = otguibase.OTStudy('Study_Crue4Stochastic')
otguibase.OTStudy.Add(Study_Crue4Stochastic)
dist_Q = ot.Gumbel(0.00179211, 1013)
dist_Q = ot.TruncatedDistribution(dist_Q, 0)
Q = otguibase.Input('Q', 1013, dist_Q, 'Débit maximal annuel (m3/s)')
dist_Ks = ot.Normal(30, 7.5)
dist_Ks = ot.TruncatedDistribution(dist_Ks, 0)
Ks = otguibase.Input('Ks', 30, dist_Ks, 'Strickler (m^(1/3)/s)')
dist_Zv = ot.Uniform(49, 51)
Zv = otguibase.Input('Zv', 50, dist_Zv, 'Côte de la rivière en aval (m)')
dist_Zm = ot.Uniform(54, 56)
Zm = otguibase.Input('Zm', 55, dist_Zm, 'Côte de la rivière en amont (m)')
Smean = otguibase.Output('Smean', '')
inputCollection = []
inputCollection.append(Q)
inputCollection.append(Ks)
inputCollection.append(Zv)
inputCollection.append(Zm)
outputCollection = []
outputCollection.append(Smean)
code='import openturns as ot\nfrom math import sqrt\n\n# 1. The function G\ndef functionCrue8vars(X) :\n    Q, Ks, Zv, Zm, Hd, Zb, L, B = X\n    Zd = Zb + Hd\n    alpha = (Zm - Zv)/L\n    H = (Q/(Ks*B*sqrt(alpha)))**(3.0/5.0)\n    Zc = H + Zv\n    S = Zc - Zd\n    return [S]\n\ndef _exec(Q, Ks, Zv, Zm):\n    # 1. Creation of the problem function\n    f8v = ot.PythonFunction(8, 1, functionCrue8vars)\n    X = [Q, Ks, Zv, Zm]\n    g = ot.ParametricFunction(f8v, [0, 1, 2, 3], X)\n    # 2. Random vector definition\n    Hd = ot.Uniform(4.,14.)\n    Zb = ot.Uniform(50.,60.)\n    L = ot.Uniform(1000.,10000.)\n    B = ot.Uniform(50.,500.)\n    inputvector = ot.ComposedDistribution([Hd, Zb, L, B])\n    inputRV = ot.RandomVector(inputvector)\n    S = ot.RandomVector(g, inputRV)\n    # 3. Sample output\n    sampleSize = 10000\n    outputSample = S.getSample(sampleSize)\n    Smean = outputSample.computeMean()[0]\n    return Smean\n'
myPhysicalModel = otguibase.PythonPhysicalModel('myPhysicalModel', inputCollection, outputCollection, code)
Study_Crue4Stochastic.add(myPhysicalModel)
tendanceCentrale_0 = otguibase.MonteCarloAnalysis('tendanceCentrale_0', myPhysicalModel)
tendanceCentrale_0.setMaximumCalls(1000)
tendanceCentrale_0.setMaximumCoefficientOfVariation(-1)
tendanceCentrale_0.setBlockSize(16)
tendanceCentrale_0.setSeed(0)
Study_Crue4Stochastic.add(tendanceCentrale_0)
