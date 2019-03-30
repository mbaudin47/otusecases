from openturns.viewer import View
import openturns as ot
from math import sqrt

# 1. The function G
def functionCrue(X) :
    Hd = 3.0
    Zb = 55.5
    L = 5.0e3
    B = 300.0
    Zd = Zb + Hd
    Q, Ks, Zv, Zm = X
    alpha = (Zm - Zv)/L
    H = (Q/(Ks*B*sqrt(alpha)))**(3.0/5.0)
    Zc = H + Zv
    S = Zc - Zd
    return [S]

# Creation of the problem function
g = ot.PythonFunction(4, 1, functionCrue) 
g = ot.MemoizeFunction(g)

# 2. Random vector definition
myParam = ot.GumbelAB(1013., 558.)
Q = ot.ParametrizedDistribution(myParam)
otLOW = ot.TruncatedDistribution.LOWER
Q = ot.TruncatedDistribution(Q, 0, otLOW)
Ks = ot.Normal(30.0, 7.5)
Ks = ot.TruncatedDistribution(Ks, 0, otLOW)
Zv = ot.Uniform(49.0, 51.0)
Zm = ot.Uniform(54.0, 56.0)

# 3. View the PDF
Q.setDescription(["Q (m3/s)"])
Ks.setDescription(["Ks (m^(1/3)/s)"])
Zv.setDescription(["Zv (m)"])
Zm.setDescription(["Zm (m)"])

View(Q.drawPDF()).show()
View(Ks.drawPDF()).show()
View(Zv.drawPDF()).show()
View(Zm.drawPDF()).show()

# 4. Create the joint distribution function, 
#    the output and the event. 
inputvector = ot.ComposedDistribution([Q, Ks, Zv, Zm])
inputRV = ot.RandomVector(inputvector)
S = ot.RandomVector(g, inputRV)
eventF = ot.Event(S, ot.GreaterOrEqual(), 0) 

# 4.bis Draw pairs
sample = inputvector.getSample(500)
myPairs = ot.Pairs(sample, "N=500", sample.getDescription(), "red", "bullet")
View(myPairs).show()

# 5. Create the Monte-Carlo algorithm
algoProb = ot.ProbabilitySimulationAlgorithm(eventF)
algoProb.setMaximumOuterSampling(1000000)
algoProb.run()

# 6. Get the results
resultAlgo = algoProb.getResult()
neval = g.getEvaluationCallsNumber()
print("Number of function calls = %d" %(neval))
pf = resultAlgo.getProbabilityEstimate()
print("Failure Probability = %e" % (pf))
level = 0.95
c95 = resultAlgo.getConfidenceLength(level)
pmin=pf-0.5*c95
pmax=pf+0.5*c95
print("%.1f %% confidence interval :[%e,%e] " % (level*100,pmin,pmax))

# 7. Plot the histogram
from openturns import VisualTest
outComputedPoints = g.getOutputHistory()
histoGraph = VisualTest.DrawHistogram(outComputedPoints,100)
histoGraph.setTitle("Histogramme de la surverse")
histoGraph.setXTitle("S (m)")
histoGraph.setYTitle("Frequence")
'''
Modifie la borne supérieure de l'axe des abscisses. 
Sinon, la surverse la plus grande est utilisée (proche de 70). 
Les valeurs de surverses plus grandes que 5 (m) sont associées 
à une densité très faible, ce qui mène à un histogramme très déformé. 
'''
bb = histoGraph.getBoundingBox()
bb_upper = bb.getUpperBound()
bb_upper[0] = 5
bb.setUpperBound(bb_upper)
histoGraph.setBoundingBox(bb)
histoGraph.setLegends([""])
View(histoGraph).show()

