import openturns as ot
import openturns.viewer as otv
import numpy as np

tmin=0. # Date minimale
tmax=12. # Date maximale
gridsize=100 # Nombre de pas de temps
mesh = ot.IntervalMesher([gridsize-1]).build(ot.Interval(tmin, tmax))

def AltiFunc(X):
    g  = 9.81
    z0 = X[0]
    v0 = X[1]
    m  = X[2]
    c  = X[3]
    zmin = X[4]
    tau=m/c
    vinf=-m*g/c
    t = np.linspace(tmin,tmax,gridsize)
    z=z0+vinf*t+tau*(v0-vinf)*(1-np.exp(-t/tau))
    z=np.maximum(z,zmin)
    altitude = [[zeta] for zeta in z]
    return altitude

inputDim = 5
outputDim = 1
alti = ot.PythonPointToFieldFunction(inputDim, mesh, outputDim, AltiFunc)


# Creation of the input distribution
distZ0 = ot.Uniform(100.0, 150.0)
distV0 = ot.Normal(55.0, 10.0)
distM = ot.Normal(80.0, 8.0)
distC = ot.Uniform(0.0, 30.0)
distZmin = ot.Dirac([0.0])
distX = ot.ComposedDistribution([distZ0, distV0, distM, distC, distZmin])

# Sample the model
size = 10
inputSample = distX.getSample(size)
outputSample = alti(inputSample)

# Draw some curves
graph = outputSample.drawMarginal(0)
graph.setTitle('chute visqueuse')
graph.setXTitle(r'$t$')
graph.setYTitle(r'$z$')
graph.setColors([ot.Drawable.ConvertFromHSV(i * (360.0/size), 1.0, 1.0) for i in range(len(graph.getDrawables()))])
ot.Show(graph)

# Compute the KL decomposition of the output
size = 500
inputSample = distX.getSample(size)
outputSample = alti(inputSample)
algo = ot.KarhunenLoeveSVDAlgorithm(outputSample, 1.0e-6)
algo.run()
KLResult = algo.getResult()
scaledModes = KLResult.getScaledModesAsProcessSample()
graph = scaledModes.drawMarginal(0)
graph.setTitle('Modes de KL, chute visqueuse')
graph.setXTitle(r'$t$')
graph.setYTitle(r'$z$')
ot.Show(graph)

# Here we have to suppress the Dirac component
distX = distX.getMarginal(range(4))
alti = ot.PointToFieldConnection(alti, ot.SymbolicFunction(["z0", "v0", "m", "c"], ["z0", "v0", "m", "c", "0.0"]))
inputSample = inputSample.getMarginal(range(4))

postProcessing = ot.KarhunenLoeveLifting(KLResult)
outputSampleChaos = KLResult.project(outputSample)

size = 20
validationInputSample = distX.getSample(size)
validationOutputSample = alti(validationInputSample)

# First, using the most basic interface
algo = ot.FunctionalChaosAlgorithm(inputSample, outputSampleChaos)
algo.run()
metaModel = ot.PointToFieldConnection(postProcessing, algo.getResult().getMetaModel())

graph = validationOutputSample.drawMarginal(0)
graph.setColors(['red'])
graph2 = metaModel(validationInputSample).drawMarginal(0)
graph2.setColors(['blue'])
graph.add(graph2)
graph.setTitle('Comparaison modele/meta-modele')
graph.setXTitle(r'$t$')
graph.setYTitle(r'$z$')
otv.View(graph)

# Second, using a more evolved interface
basis = ot.OrthogonalProductPolynomialFactory([ot.StandardDistributionPolynomialFactory(distX.getMarginal(i)) for i in range(distX.getDimension())])
adaptiveStrategy = ot.FixedStrategy(basis, ot.EnumerateFunction(distX.getDimension()).getStrataCumulatedCardinal(6))
projectionStrategy = ot.LeastSquaresStrategy(ot.LeastSquaresMetaModelSelectionFactory(ot.LARS(), ot.CorrectedLeaveOneOut()))
algo = ot.FunctionalChaosAlgorithm(inputSample, outputSampleChaos, distX, adaptiveStrategy, projectionStrategy)
algo.run()
metaModel = ot.PointToFieldConnection(postProcessing, algo.getResult().getMetaModel())

graph = validationOutputSample.drawMarginal(0)
graph.setColors(['red'])
graph2 = metaModel(validationInputSample).drawMarginal(0)
graph2.setColors(['blue'])
graph.add(graph2)
graph.setTitle('Comparaison modele/meta-modele')
graph.setXTitle(r'$t$')
graph.setYTitle(r'$z$')
otv.View(graph)
