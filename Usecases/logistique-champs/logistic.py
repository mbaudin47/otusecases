import openturns as ot
from openturns.viewer import View

# Here the model is created using the simplified Python interface for FieldToPointFunction

tmin=1790. # Date minimale
tmax=2001. # Date maximale
gridsize=100 # Nombre de pas de temps
modeleName = "Modèle Logistique" # Nom du modèle
parameterIndexName = "Temps (années)" # Nom du paramètre d'indexation
fieldName = "Population (millions)" # Nom du champ
alphaInf = 0.05 # Valeur du niveau alpha pour le quantile inférieur
# Définit l'intervalle de temps pour la simulation
monhorizon = ot.Interval(tmin, tmax)
# Définit la grille temporelle régulière
mesh = ot.IntervalMesher([gridsize-1]).build(monhorizon)
graph1 = mesh.draw()
graph1.setTitle(modeleName)
graph1.setXTitle(parameterIndexName)
View(graph1)

from numpy import exp
from numpy import array

def logisticSolution(X):
    # Récupère les noeuds du maillage
    v = mesh.getVertices()
    # Convertit en tableau
    t = array(v)
    # Convertit le tableau 2D en tableau 1D
    t = t.flatten()
    # Récupère la date initiale
    t0 = t[0]
    # Calcule la trajectoire
    y0,a,b = X
    y=a*y0/(b*y0+(a-b*y0)*exp(-a*(t-t0)))
    y = y/1.e6
    # Créée la liste des altitudes à partir du array numpy
    matrajectoire = [[zeta] for zeta in y]
    return matrajectoire

inputDim = 3 # Nombre de variables en entrée (i.e. taille du vecteur aléatoire)
outputDim = 1 # Nombre de champs en sortie (ici, une seule trajectoire)
maFonctionChamp = ot.PythonPointToFieldFunction(inputDim, mesh, outputDim, logisticSolution)

# Teste une évaluation
y0=3.9e6 # Population initiale
a=0.03134
b=1.5887e-10
X = [y0,a,b]
unetrajectoire = maFonctionChamp(X)
unetrajectoire.setDescription([fieldName])
#print("Une trajectoire:")
#print(unetrajectoire)

# Creation of the input distribution
distY0 = ot.Normal(y0, 0.1 * y0)
distA  = ot.Normal(a, 0.3 * a)
distB  = ot.Normal(b, 0.3 * b)
distX = ot.ComposedDistribution([distY0, distA, distB])

# Sample the model
size = 10
inputSample = distX.getSample(size)
outputSample = maFonctionChamp(inputSample)
# outputSample is a ProcessSample

# Draw some trajectories
graph = outputSample.drawMarginal(0)
graph.setTitle(modeleName)
graph.setXTitle(parameterIndexName)
graph.setYTitle(fieldName)
myTrajectories =  [ot.Drawable.ConvertFromHSV(i * (360.0/size), 1.0, 1.0) for i in range(len(graph.getDrawables()))]
graph.setColors(myTrajectories)
View(graph)

# Dessine la trajectoire moyenne
size = 100
inputSample = distX.getSample(size)
outputSample = maFonctionChamp(inputSample)
champsMoyen = outputSample.computeMean()
graphMoy = champsMoyen.draw()
graphMoy.setColors(["green"])
quantileSup = outputSample.computeQuantilePerComponent(1-alphaInf)
quantileInf = outputSample.computeQuantilePerComponent(alphaInf)
graphSup = quantileSup.draw()
graphSup.setColors(["blue"])
graphInf = quantileInf.draw()
graphInf.setColors(["red"])
graphMoy.add(graphSup)
graphMoy.add(graphInf)
graphMoy.setTitle('Trajectoire')
graphMoy.setXTitle(parameterIndexName)
graphMoy.setYTitle(fieldName)
graphMoy.setLegends(["Moyenne","Quantile %.2f%%" % (100-alphaInf*100),"Quantile %.2f%%" % (alphaInf*100)])
graphMoy.setLegendPosition("topleft")
View(graphMoy)

# Compute the KL decomposition of the output
size = 100
inputSample = distX.getSample(size)
outputSample = maFonctionChamp(inputSample)
threshold = 1.e-5 # Seuil pour la troncature des valeurs propres
algo = ot.KarhunenLoeveSVDAlgorithm(outputSample, threshold)
algo.run()
KLResult = algo.getResult()
scaledModes = KLResult.getScaledModesAsProcessSample()
nbModes = scaledModes.getSize()

# Plot the KL decomposition
graph = scaledModes.drawMarginal(0)
graph.setTitle('%s, seuil=%.2e, %d modes de KL' % (modeleName,threshold,nbModes))
graph.setXTitle(parameterIndexName)
graph.setYTitle(fieldName)
modesStr = ["Mode "+str(i) for i in range(nbModes)]
graph.setLegends(modesStr)
graph.setLegendPosition("topleft")
View(graph)

# Valeurs propres selectionnees
eigenValues = KLResult.getEigenValues()
print("Nombre de modes = %d" % (eigenValues.getSize()))
print("eigenvalues=")
print(eigenValues)

# Graphe des valeurs propres selectionnees
graph = ot.Graph()
mesVP = ot.Cloud(range(eigenValues.getSize()), eigenValues)
graph.add(mesVP)
graph.setGrid(True)
graph.setTitle("%s, seuil=%.2e, %d modes de KL" % (modeleName,threshold,nbModes))
graph.setXTitle("Index")
graph.setYTitle("Eigenvalue")
graph.setAxes(True)
View(graph)

# Graphe des valeurs propres selectionnees
from numpy import cumsum
evcs = cumsum(array(eigenValues))
evcs = evcs/evcs[-1]
graph = ot.Graph()
mesVP = ot.Cloud(range(len(evcs)), evcs)
graph.add(mesVP)
graph.setGrid(True)
graph.setTitle("%s, seuil=%.2e, %d modes de KL" % (modeleName,threshold,nbModes))
graph.setXTitle("Index")
graph.setYTitle("Cumulated Eigenvalue Sum")
graph.setAxes(True)
View(graph)

# Mean function
# Champ moyen: field
mean = outputSample.computeMean()

# Tendance moyenne approximee lineaire entre les points
# appel via une methode P1
times = mesh.getVertices()
# Créée une liste des dates
locations = [t[0] for t in times]
values = mean.getValues()
myEvaluation = ot.PiecewiseLinearEvaluation(locations, values)
graph = myEvaluation.draw(tmin,tmax,100)
graph.setTitle("Mean trajectory")
graph.setXTitle(parameterIndexName)
graph.setYTitle(fieldName)
View(graph)

# Correlation function
cov = KLResult.getCovarianceModel()
cov.setName("chute visqueuse")
asCorrelation = True
isStationary = False
graph = cov.draw(0, 0, tmin, tmax, 128, isStationary, asCorrelation)
View(graph)

# Sample of coefficients Xi
sampleKsi = KLResult.project(outputSample)

# Chaque marginale est reconstruite par noyau gaussien
# False, 0, False: pas de binning (non aggergation des donnees dnas des segments), nbre de bins, pas d'effet de bord
nbmodes = sampleKsi.getDimension()
xi_marges = [ot.KernelSmoothing(ot.Normal(), False, 0, False).build(sampleKsi.getMarginal(i)) for i in range(nbmodes)]

# graphes des pdf marginales
for i in range(len(xi_marges)):
    monHisto = ot.HistogramFactory().build(sampleKsi.getMarginal(i)).drawPDF()
    monHisto.setColors(["blue"])
    graph.add(monHisto)
    graph = xi_marges[i].drawPDF()
    graph.setTitle(r"Distribution of $\xi$" + str(i))
    graph.setXTitle(r"$\xi$" + str(i))
    graph.setYTitle("PDF")
    graph.add(monHisto)
    graph.setLegends(["KS","Histogram"])
    View(graph)

# graphe des pairs dans l'espace des rangs
pairs = ot.Pairs(sampleKsi.rank())
pairs.setPointStyle("bullet")
View(pairs)

