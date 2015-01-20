from evaluation import *
from evaluationSparsifiers import *
from evaluationProperties import *

#Returns a list of target edge ratios
def getEdgeRatios():
	#return [1.0]
	return [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
	#return [0.01, 0.1, 0.3, 0.6, 0.8, 0.9]

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
        S_ForestFire("1", 0.15, 5)
        ]

#Returns a list of graphs
def getGraphs():
    return [
        #GraphDescription("./input/Yale4.graphml", Format.GraphML, "Yale4"),
        #GraphDescription("./input/Virginia63.graphml", Format.GraphML, "Virginia63"),
        #GraphDescription("./input/Tennessee95.graphml", Format.GraphML, "Tennessee95"),
        GraphDescription("./input/us-aviation-t100-2013.graphml", Format.GraphML, "USAviation"),
        GraphDescription("./input/Caltech36.graphml", Format.GraphML, "Caltech36"),
        GraphDescription("./input/kitEmail.graphml", Format.GraphML, "KitEmail"),
        GraphDescription("./input/LFR-1000.graph", Format.METIS, "LFR-1000"),
        GraphDescription("./input/PGPgiantcompo.graph", Format.METIS, "PGP"),
        GraphDescription("./input/karate.graph", Format.METIS, "Karate"),
        GraphDescription("./input/bter-graph-cc0.7-small.graph", Format.EdgeListCommaOne, "BTER"),
        GraphDescription("./input/ErdosRenyi.graph", Format.METIS, "ErdosRenyi"),
        GraphDescription("./input/jazz.graph", Format.METIS, "Jazz")
        ]

#Returns a list of property calculation descriptions
def getProperties():
    return [
        P_General(),
        P_Ratios(),
        #P_Community(),
        #P_Diameter(),
        #P_DegreeDistribution(),
        #P_Centrality(),
        #P_Components()
		P_KolmogorowSmirnow()
        ]
