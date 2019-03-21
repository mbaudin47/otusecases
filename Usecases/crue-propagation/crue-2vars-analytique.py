#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import openturns as ot
import otguibase

myStudy = otguibase.OTStudy('Study_Crue2Vars')

dist_Q = ot.TruncatedDistribution(ot.Gumbel(1./558., 1013.), 0, ot.TruncatedDistribution.LOWER)
dist_Ks = ot.TruncatedDistribution(ot.Normal(30.0, 7.5), 0, ot.TruncatedDistribution.LOWER)

Q = otguibase.Input('Q', 1000., dist_Q, 'Débit maximal annuel (m3/s)')
Ks = otguibase.Input('Ks', 30., dist_Ks, 'Strickler (m^(1/3)/s)')
H = otguibase.Output('H', 'Hauteur d\'eau (m)')

inputs = [Q, Ks]
outputs = [H]
formulas = ['(Q/(Ks*300.*sqrt(0.001)))^(3./5.)']
myPhysicalModel = otguibase.SymbolicPhysicalModel('myPhysicalModel', inputs, outputs, formulas)
myStudy.add(myPhysicalModel)

## script
script = myStudy.getPythonScript()
print(script)
exec(script)
