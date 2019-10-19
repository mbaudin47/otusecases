from openturns.viewer import View
import numpy as np
import openturns as ot

# 1. The function G
def modeleChaboche(X):
    strain,R,C,gamma = X
    stress = R + C*(1-np.exp(-gamma*strain))
    return [stress]

# Creation of the problem function
f = ot.PythonFunction(4, 1, modeleChaboche) 

# 2. Random vector definition
Strain = ot.Uniform(0,0.07)
unknownR = 750e6
unknownC = 2750e6
unknownGamma = 10
R = ot.Dirac(unknownR)
C = ot.Dirac(unknownC)
Gamma = ot.Dirac(unknownGamma)

# 3. View the PDF
Strain.setDescription(["Strain"])
R.setDescription(["R"])
C.setDescription(["C"])
Gamma.setDescription(["Gamma"])

# 4. Create the joint distribution function
inputRandomVector = ot.ComposedDistribution([Strain, R, C, Gamma])

# 5. Create the Monte-Carlo algorithm
sampleSize = 100
inputSample = inputRandomVector.getSample(sampleSize)
#print(inputSample)
outputSigma = f(inputSample)
#print(outputSigma)

# 7. Plot the histogram
histoGraph = ot.HistogramFactory().build(outputSigma/1.e6).drawPDF()
histoGraph.setTitle("Histogramme de la contrainte")
histoGraph.setXTitle("Stress (MPa)")
histoGraph.setYTitle("Frequence")
#histoGraph.setBoundingBox([-1,7,0,0.60])
histoGraph.setLegends([""])
View(histoGraph)

# Generate observation noise
sigmaObservationNoiseSigma = 40.e6 # (Pa)
noiseSigma = ot.Normal(0.,sigmaObservationNoiseSigma)
sampleNoiseH = noiseSigma.getSample(sampleSize)
observedSigma = outputSigma + sampleNoiseH

# Create and save sample
observedSample = ot.Sample(sampleSize,2)
observedSample.setDescription(["Strain","Stress"])
observedSample[:,0] = inputSample[:,0]
observedSample[:,1] = observedSigma[:]

observedSample.exportToCSVFile("chaboche-observations.csv")

