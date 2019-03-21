#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import openturns as ot
import otguibase
from inspect import getframeinfo, currentframe
from os.path import join, dirname, abspath
import otmorris

myStudy = otguibase.OTStudy('Study_Morris')

## Model
dist_X1 = ot.Uniform(0., 1.)

Y = otguibase.Output('Y','Y')

#
# Parameters of the Morris function
b0 = ot.DistFunc.rNormal()
alpha = ot.DistFunc.rNormal(10)
beta = ot.DistFunc.rNormal(6*14)
gamma = ot.DistFunc.rNormal(20*14)
myMF = otmorris.MorrisFunction(alpha, beta, gamma, b0)
myMorrisFunction = ot.NumericalMathFunction( myMF )

#
# Python function code
code = '''
def _exec(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20):
    X = [X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20]
    Y = myMorrisFunction(X)[0]
    return Y
'''

# Define the physical model
d = 20
inputs = []
for i in range(d):
    namei = "X"+str(i+1)
    desci = ""
    Xi = otguibase.Input(namei, 0.5, dist_X1, desci)
    inputs.append(Xi)
model = otguibase.PythonPhysicalModel('myPhysicalModel', inputs, [Y], code)
myStudy.add(model)

## script
script = myStudy.getPythonScript()
print(script)
exec(script)
