# -*- coding: utf-8 -*-
# Copyright (C) 2018 - 2019 - Michael Baudin

import pylab as pl
import openturns as ot
import numpy as np

'''
Loi de Chaboche pour l’endommagement par fatigue. 
C'est une loi d'écrouissage qui caractérise l'évolution de 
la contrainte en fonction de la déformation. 
Remerciements : A.Dumas, T.Yalamas (Phiméca)
'''

# Donnees reelles
observedSample = ot.Sample_ImportFromCSVFile("chaboche-observations-v1.csv")
nbobs = observedSample.getSize()
strainObservations = observedSample[:,0]
stressObservations = observedSample[:,1]

# Parametres
R = 750e6 # Exact : 750e6
C = 2700e6 # Exact : 2750e6
Gamma = 8. # Exact : 10
theta0 = [R,C,Gamma]

# define the problem bounds
boundsMin = [600e6,2000e6,7.]
boundsMax = [800e6,3000e6,12.]

# Modèle
def modeleChaboche(X):
    strain,R,C,Gamma = X
    stress = R + C*(1-np.exp(-Gamma*strain))
    return [stress]

# Fonction de calibration
def calibrationG(theta):
    R,C,Gamma = theta
    predictedStress = ot.Sample(nbobs,1)
    for i in range(nbobs):
        X = [strainObservations[i,0],R,C,Gamma]
        predictedStress[i,0] = modeleChaboche(X)[0]
    return predictedStress.asPoint()

calibrationFunc = ot.PythonFunction(3, nbobs, calibrationG)

def costFunction(theta):
    # Fonction coût pour l'assimilation de données 
    # Calcule les résidus
    stressPredicted = calibrationFunc(theta)
    r = stressObservations.asPoint() - stressPredicted
    # Fait la somme des carrés
    sumOfSquares = r.normSquare()**2
    return [sumOfSquares]

dimCalage = len(theta0)
objective = ot.PythonFunction(dimCalage, 1, costFunction)

bounds = ot.Interval(boundsMin,boundsMax)
# define the problem
problem = ot.OptimizationProblem(objective)
problem.setMinimization(True)
problem.setBounds(bounds)

labelsTheta = ("R","C","Gamma")

# Résout le NLLS
maximumIteration = 1000
algo = ot.Cobyla()
algo.setProblem(problem)
algo.setMaximumIterationNumber(maximumIteration)
algo.setStartingPoint(theta0)
# Lance l'optimisation
algo.run()
# retrieve results
result = algo.getResult()
thetaStar = result.getOptimalPoint()
# Get the number of function evaluations
nfeval = objective.getEvaluationCallsNumber()

# Pour vérifier que tout est OK 
print("Coût initial:")
print(costFunction(theta0))

# Affiche le coût avant optimisation
dimCalage = len(theta0)
costInitial = costFunction(theta0)
print("Theta0=")
for i in range(dimCalage):
    print("\t%s 0[%d]=%.3e" % (labelsTheta[i],i,theta0[i]))
print("Cout initial=%.3e" % (costInitial[0]))

# Affiche l'optimum
print('Theta*=')
for i in range(dimCalage):
    print("\t%s *[%d]=%.3e" % (labelsTheta[i],i,thetaStar[i]))
# Coût optimal
costFinal = costFunction(thetaStar)
print("Cout final=%.3e" % (costFinal[0]))
print("Nombre d'évaluations=%d" % (nfeval))

# Plot the fit before/after
pl.xlabel("Strain")
pl.ylabel("Stress (Pa)")
pl.title("Calage NLLS")
pl.plot(strainObservations,stressObservations,"bo",label="Data")
y=calibrationFunc(theta0)
pl.plot(strainObservations,y,"ro",label="Before Calibration")
y=calibrationFunc(thetaStar)
pl.plot(strainObservations,y,"go",label="After Calibration")
pl.legend()
