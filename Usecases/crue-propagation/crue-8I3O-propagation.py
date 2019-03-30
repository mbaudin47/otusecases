"""
Propagation des incertitudes 
à travers un modèle de la hauteur d'eau d'une rivière. 

Reference:
Derivative-based global sensitivity measures: general links
with Sobol’ indices and numerical tests
M. Lamboni, B. Iooss, A.-L. Popelin, F. Gamboa

Revue sur l’analyse de sensibilite globale de modeles
numériques, Bertrand Iooss, 7 Dec 2010
"""

import openturns as ot
from math import sqrt, expm1
from openturns.viewer import View

# 1. Define the G function
def functionCrue(X) :
    Q, Ks, Zv, Zm, Hd, Zb, L, B = X
    alpha = (Zm - Zv)/L
    H = (Q/(Ks*B*sqrt(alpha)))**(3.0/5.0)
    Zc = H + Zv
    Zd = Zb + Hd
    S = Zc - Zd
    if (S<0):
        CS = 0.2+0.8-expm1(-1000/S**4)
    else:
        CS = 1
    if (Hd<8):
        CH = 8./20.
    else:
        CH = Hd/20.
    C=CS+CH
    return [H,S,C]

myFunction = ot.PythonFunction(8, 3, functionCrue) 

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

# 3. Create the joint distribution
inputDistribution = ot.ComposedDistribution((Q,Ks,Zv,Zm,Hd,Zb,L,B))
inputRandomVector = ot.RandomVector(inputDistribution)
outputRandomVector = ot.RandomVector(myFunction, inputRandomVector)

# 4. Simple Monte-Carlo sampling
samplesize=500
outputSample=outputRandomVector.getSample(samplesize)

# 6. Plot the histogram
from openturns import VisualTest
histoGraph = VisualTest.DrawHistogram(outputSample[:,0])
histoGraph.setTitle("Histogramme de la hauteur")
histoGraph.setXTitle("H (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setLegends([""])
View(histoGraph).show()

histoGraph = VisualTest.DrawHistogram(outputSample[:,1])
histoGraph.setTitle("Histogramme de la surverse")
histoGraph.setXTitle("S (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setLegends([""])
View(histoGraph).show()
