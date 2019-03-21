#!/usr/bin/python
#coding=utf8

"""
Polynomial chaos expansions on a nonlinear oscillator.

This script is compatible with OpenTurns = 1.12

Authors:
    Géraud Blatman (EDF R&D)
    Chu Mai (EDF R&D)
    Michaël Baudin
"""

# Import useful modules
import numpy as np
import openturns as ot
# set seed for random number generator to assure the reproducibility of the analysis
ot.RandomGenerator.SetSeed(1)

# Nonlinear oscillator function
#
dim = 7
#
def NL_oscillator(x):
    mp, ms, kp, ks, xip, xis, S0 = x
    omega_p = np.sqrt(kp/mp)
    omega_s = np.sqrt(ks/ms)
    gamma = ms/mp
    omega_a = 0.5*(omega_p+omega_s)
    xi_a = 0.5*(xip+xis)
    theta = 1./omega_a*(omega_p-omega_s)
    F = \
        3*ks*np.sqrt(np.pi*S0/(4.*xis*omega_s**3)*
        xi_a*xis/(xip*xis*(4.*xi_a**2+theta**2)+gamma*xi_a**2)*
        (xip*omega_p**3+xis*omega_s**3)*omega_p/(4.*xi_a*omega_a**4))
    return [F]
#
f = ot.PythonFunction(dim,1,NL_oscillator)


## Define the input probabilistic model
mean_list = [1.5, 0.01, 1., 0.01, 0.05, 0.02, 100.]
cov_list = [0.1, 0.1, 0.2, 0.2, 0.4, 0.5, 0.1]
myCollection = ot.DistributionCollection(dim)
for i, (mu, cov) in enumerate(zip(mean_list, cov_list)):
    myParam = ot.LogNormalMuSigma(mu, mu*cov, 0.)
    myCollection[i] = ot.ParametrizedDistribution(myParam)

distribution = ot.ComposedDistribution(myCollection)
rv = ot.RandomVector(distribution)

## Construct a polyomial chaos approximation

n = 1000
p = 4

# Polynomial basis
polyColl = ot.PolynomialFamilyCollection(dim)
for i in range(dim):
    polyColl[i] = \
        ot.StandardDistributionPolynomialFactory(distribution.getMarginal(i))
enumerateFunction = ot.LinearEnumerateFunction(dim)
basis = ot.OrthogonalProductPolynomialFactory(polyColl, enumerateFunction)


## Get input and output samples
#sampleX = LowDiscrepancyExperiment(SobolSequence(), distribution, n).generate()
sampleX = ot.LHSExperiment(distribution, n).generate()
sampleY = f(sampleX)
weights = ot.Point(n ,1) # vector of unit weights: [1,1,...,1]
basisSequenceFactory = ot.LARS()
fittingAlgorithm = ot.CorrectedLeaveOneOut()
approximationAlgorithm = ot.LeastSquaresMetaModelSelectionFactory(basisSequenceFactory, fittingAlgorithm)
# Construct the chaos approximation
P = enumerateFunction.getStrataCumulatedCardinal(p) # basis size
fixedStrategy = ot.FixedStrategy(basis, P)
evalStrategy = ot.LeastSquaresStrategy(sampleX, weights, sampleY,approximationAlgorithm)
algo = ot.FunctionalChaosAlgorithm(sampleX, weights, sampleY, \
                                distribution, fixedStrategy, evalStrategy)
algo.run()

# Get the chaos result
result = algo.getResult()
metamodel = result.getMetaModel()

## Check the chaos accuracy using a validation basis

# Get validation input and output samples
n_valid = n*3
sampleX_valid = ot.LHSExperiment(distribution, n_valid).generate()
sampleY_valid = f(sampleX_valid)

# Chaos predictions
sampleY_pred = metamodel(sampleX)
sampleY_valid_pred = metamodel(sampleX_valid)

# Compute the error estimates R2 and Q2
residuals = sampleY -sampleY_pred
residuals_pred = sampleY_valid-sampleY_valid_pred
R2 = 1. - residuals.computeVariance()[0] / sampleY.computeVariance()[0]
Q2 = 1. - residuals_pred.computeVariance()[0] / sampleY_valid.computeVariance()[0]

# Plot the actual vs the predicted outputs
title = "Approximation accuracy  -  R2=%.2f  -  Q2=%.2f" % (R2, Q2)
title = "Approximation accuracy  -  R2=%.2f  -  Q2=%.2f" % (R2, Q2)
valid = ot.MetaModelValidation(sampleX_valid, sampleY_valid, metamodel)
graph = valid.drawValidation()
graph.setTitle(title)
ot.viewer.View(graph).show()

## Compute the mean and the standard deviation of the model output

# Monte Carlo (sampling) approach
bigSampleX = rv.getSample(10000)
bigSampleY = f(bigSampleX)
mean_mc = bigSampleY.computeMean()[0]
std_mc = bigSampleY.computeStandardDeviationPerComponent()[0]
print("Monte Carlo - Mean: %.2f, StD: %.2f" % (mean_mc, std_mc))

# Chaos-based estimates
chaosRV = ot.FunctionalChaosRandomVector(result)
mean_mc = chaosRV.getMean()[0]
std_mc = np.sqrt(chaosRV.getCovariance()[0,0])
print("Chaos - Mean: %.2f, StD: %.2f" % (mean_mc, std_mc))


## Compute the Sobol indices
sobol_indices = np.array([chaosRV.getSobolIndex(i) for i in range(dim)])
sobol_total_indices = np.array([chaosRV.getSobolTotalIndex(i) for i in range(dim)])
print("\nSobol indices: \n", sobol_indices)
print("\nTotal sobol indices: \n", sobol_total_indices)

