from openturns.viewer import View
from numpy import inf
import openturns as ot
from math import sqrt

# 1. The function G
def functionCrue(X) :
    L = 5.0e3
    B = 300.0
    Q, K_s, Z_v, Z_m = X
    alpha = (Z_m - Z_v)/L
    H = (Q/(K_s*B*sqrt(alpha)))**(3.0/5.0)
    return [H]

# Creation of the problem function
f = ot.PythonFunction(4, 1, functionCrue) 
f.enableHistory()

# 2. Random vector definition
Q = ot.Gumbel(1./558., 1013.)
print(Q)

'''
Q = ot.Gumbel()
Q.setParameter(ot.GumbelAB()([1013., 558.]))
print(Q)
'''
Q = ot.TruncatedDistribution(Q, 0, inf)
unknownKs = 30.0
unknownZv = 50.0
unknownZm = 55.0
K_s = ot.Dirac(unknownKs)
Z_v = ot.Dirac(unknownZv)
Z_m = ot.Dirac(unknownZm)

# 3. View the PDF
Q.setDescription(["Q (m3/s)"])
K_s.setDescription(["Ks (m^(1/3)/s)"])
Z_v.setDescription(["Zv (m)"])
Z_m.setDescription(["Zm (m)"])

# 4. Create the joint distribution function
inputRandomVector = ot.ComposedDistribution([Q, K_s, Z_v, Z_m])

# 5. Create the Monte-Carlo algorithm
sampleSize = 100
inputSample = inputRandomVector.getSample(sampleSize)
#print(inputSample)
outputH = f(inputSample)
#print(outputH)

# 7. Plot the histogram
from openturns import VisualTest
histoGraph = VisualTest.DrawHistogram(outputH,20)
histoGraph.setTitle("Histogramme de la surverse")
histoGraph.setXTitle("S (m)")
histoGraph.setYTitle("Frequence")
histoGraph.setBoundingBox([-1,7,0,0.60])
histoGraph.setLegends([""])
View(histoGraph).show()

# Generate observation noise
sigmaObservationNoiseH = 0.1 # (m)
noiseH = ot.Normal(0.,sigmaObservationNoiseH)
sampleNoiseH = noiseH.getSample(sampleSize)
observedH = outputH + sampleNoiseH

# Create and save sample
observedSample = ot.Sample(sampleSize,2)
observedSample.setDescription(["Q (m3/s)","H (m)"])
observedSample[:,0] = inputSample[:,0]
observedSample[:,1] = observedH[:]

observedSample.exportToCSVFile("crue-debit-hauteur.csv")
