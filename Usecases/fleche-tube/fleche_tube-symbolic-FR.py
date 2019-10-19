# -*- coding: utf-8 -*-
"""
Faire une étude en tendance centrale sur la formule de 
la flèche d'un tube.
"""

import openturns as ot
from openturns.viewer import View

# 1. Define the function
inputsvars=["F","L","a","De","di","E"]
formula = ["-F*a^2*(L-a)^2/(3*E*L*pi_*(De^4-di^4)/32)"]
g = ot.SymbolicFunction(inputsvars,formula)

# 2. Define the probabilistic model
XF=ot.Normal(1,0.1)
XL=ot.Normal(1.5,0.01)
Xa=ot.Uniform(0.7,1.2)
XD=ot.Triangular(0.75,0.8,0.85)
Xd=ot.Triangular(0.09,0.1,0.11)
XE=ot.Normal(200000,2000)
XF.setDescription(["Sollicitation"])
XL.setDescription(["Longueur"])
Xa.setDescription(["Appui"])
XD.setDescription(["Diametre externe"])
Xd.setDescription(["Diametre interne"])
XE.setDescription(["Module d'Young"])

# 3. Create the joint distribution function, 
#    the output and the event. 
inputDistribution = ot.ComposedDistribution([XF,XL,Xa,XD,Xd,XE])

# Create the input random vector
inputRandomVector = ot.RandomVector(inputDistribution)

outputRV = ot.CompositeRandomVector(g, inputRandomVector)

# 7. Plot the histogram
outputSample = outputRV.getSample(500)
histoGraph = ot.HistogramFactory().build(outputSample).drawPDF()
samplesize=outputSample.getSize()
histoGraph.setTitle("Histogramme de la flèche - N=%d" % (samplesize))
histoGraph.setXTitle("Flèche")
histoGraph.setYTitle("Fréquence")
histoGraph.setLegends([""])
View(histoGraph)

