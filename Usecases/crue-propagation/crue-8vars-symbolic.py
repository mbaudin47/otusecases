#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import openturns as ot
import otguibase

myStudy = otguibase.OTStudy('Study_Crue8vars')

## Model
dist_Q = ot.TruncatedDistribution(ot.Gumbel(1./558., 1013.), 0, ot.TruncatedDistribution.LOWER)
dist_Ks = ot.TruncatedDistribution(ot.Normal(30.0, 7.5), 0, ot.TruncatedDistribution.LOWER)
dist_Zv = ot.Uniform(49.0, 51.0)
dist_Zm = ot.Uniform(54.0, 56.0)
dist_Hd = ot.Uniform(7, 9)
dist_Zb = ot.Triangular(55, 55.5, 56)
dist_L = ot.Triangular(4990, 5000, 5010)
dist_B = ot.Triangular(295, 300, 305)

Q = otguibase.Input('Q', 1013., dist_Q, 'Débit maximal annuel (m3/s)')
Ks = otguibase.Input('Ks', 30., dist_Ks, 'Strickler (m^(1/3)/s)')
Zv = otguibase.Input('Zv', 50., dist_Zv, 'Côte de la rivière en aval (m)')
Zm = otguibase.Input('Zm', 55., dist_Zm, 'Côte de la rivière en amont (m)')
Hd = otguibase.Input('Hd', 8, dist_Hd, 'Hauteur de la digue (m)')
Zb = otguibase.Input('Zb', 55.5, dist_Zb, 'Hauteur de la berge (m)')
L = otguibase.Input('L', 5000, dist_L, 'Longueur de la rivière (m)')
B = otguibase.Input('B', 300, dist_B, 'Largeur de la rivière (m)')

S = otguibase.Output('S','Surverse (m)')

inputs = [Q,Ks,Zv,Zm,Hd,Zb,L,B]
outputs = [S]
formulas = ['(Q/(Ks*B*sqrt((Zm-Zv)/L)))^(3.0/5.0)+Zv-Zb-Hd']
myPhysicalModel = otguibase.SymbolicPhysicalModel('myPhysicalModel', inputs, outputs, formulas)
myStudy.add(myPhysicalModel)

## script
script = myStudy.getPythonScript()
print(script)
exec(script)
