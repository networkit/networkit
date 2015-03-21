from evaluation import *
from evaluationSparsifiers import *
from evaluationProperties import *

#Returns a list of target edge ratios
def getEdgeRatios():
	#return [1.0]
	return [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	#return [0.01, 0.1, 0.3, 0.6, 0.8, 0.9]
	#return [0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]

#Returns a list of backbone algorithm descriptions
def getAlgorithms():
    return [
	S_Original(),
	S_SimmelianBackboneNonParametric(),
	S_SimmelianBackboneParametric(10),
    S_LocalSimilarity(),
	S_SimmelianMultiscale(),
	S_Multiscale(),
	S_Random(""),
	S_LocalDegree(),
	S_ForestFire("", 0.7, 5),
	S_DegreeMultiscale("lambda dx,dy: max(dx, dy)", "max"),
	S_DegreeMultiscale("lambda dx,dy: min(dx, dy)", "min"),
	S_DegreeMultiscale("lambda dx,dy: (dx + dy / 2.0)", "avg")
    ]

#Returns a list of graphs
def getGraphs():
    return [
	#thesis set
	#GraphDescription("./input/Yale4.graphml", "Yale4", Format.GraphML),
	#GraphDescription("./input/Virginia63.graphml", "Virginia63", Format.GraphML),
	#GraphDescription("./input/Tennessee95.graphml", "Tennessee95", Format.GraphML),
	GraphDescription("./input/us-aviation-t100-2013.graphml", "USAviation", Format.GraphML),
	GraphDescription("./input/Caltech36.graphml", "Caltech36", Format.GraphML),
	GraphDescription("./input/kitEmail.graphml", "KitEmail", Format.GraphML),
	GraphDescription("./input/LFR-1000.graph", "LFR-1000", Format.METIS),
	GraphDescription("./input/PGPgiantcompo.graph", "PGP", Format.METIS),
	GraphDescription("./input/karate.graph", "Karate", Format.METIS),
	GraphDescription("./input/bter-graph-cc0.7-small.graph", "BTER", Format.EdgeListCommaOne),
	GraphDescription("./input/ErdosRenyi.graph", "ErdosRenyi", Format.METIS),
	GraphDescription("./input/jazz.graph", "Jazz", Format.METIS),
	#sampling from large graphs
	#GraphDescription("./input/cit-HepTh.edgelist-t0.graph", "HepTh", Format.EdgeList, continuous=False, separator='\t', firstNode=0),
	#GraphDescription("./input/cit-HepPh.edgelist-t0.graph", "HepPh", Format.EdgeList, continuous=False, separator='\t', firstNode=0),
	#GraphDescription("./input/soc-Epinions1.edgelist-t0.graph", "Epinions", Format.EdgeListTabZero),
	#GraphDescription("./input/as20000102.txt", "AS", Format.EdgeList, continuous=False, separator='\t', firstNode=1),
	#some more fb graphs
	#GraphDescription("./input/fb-FSU53.edgelist", "fb-FSU53", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-Indiana69.edgelist", "fb-Indiana69", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-Michigan23.edgelist", "fb-Michigan23", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-MSU24.edgelist", "fb-MSU24", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-Penn94.edgelist", "fb-Penn94", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-Texas80.edgelist", "fb-Texas80", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-Texas84.edgelist", "fb-Texas84", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-UF21.edgelist", "fb-UF21", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-UGA50.edgelist", "fb-UGA50", Format.EdgeListTabZero),
	#GraphDescription("./input/fb-UIllinois20.edgelist", "fb-UIllinois20", Format.EdgeListTabZero),
	#others
	#GraphDescription("./input/mouse_brain_1.graphml", "mouse_brain_1", Format.GraphML),
	#GraphDescription("./input/test.fiber.small.graphml", "test.fiber.small", Format.GraphML),
	#web
	#GraphDescription("./input/cnr-2000.metis.graph", "cnr-2000", Format.METIS),
	#big..
	#GraphDescription("./input/uk-2002.metis.graph", "uk-2002", Format.METIS)
	#GraphDescription("./input/eu-2005.metis.graph", "eu-2005", Format.METIS),
	#GraphDescription("./input/in-2004.metis.graph", "in-2004", Format.METIS),
	#GraphDescription("./input/con-fiber_big.metis.graph", "fiber_big", Format.METIS),
]

#Returns a list of property calculation descriptions
def getProperties():
    return [
	P_General(),
	P_Ratios(),
	P_DegreeDistribution(),
	P_ClusteringCoefficients(),
	P_Diameter(),
	P_ConnectedComponents(),
	P_Community(),
	P_PageRank(),
	P_Betweenness()
	]
