from evaluation import *
from evaluationSparsifiers import *
from evaluationProperties import *

#Returns a list of target edge ratios
def getEdgeRatios():
	#return [1.0]
	return [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
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
        S_ForestFire("", 0.5, 5)
        ]

#Returns a list of graphs
def getGraphs():
    return [
        #GraphDescription("./input/Yale4.graphml", Format.GraphML, "Yale4"),
        #GraphDescription("./input/Virginia63.graphml", Format.GraphML, "Virginia63"),
        #GraphDescription("./input/Tennessee95.graphml", Format.GraphML, "Tennessee95"),
        #GraphDescription("./input/us-aviation-t100-2013.graphml", "USAviation", Format.GraphML),
        #GraphDescription("./input/Caltech36.graphml", "Caltech36", Format.GraphML),
        #GraphDescription("./input/kitEmail.graphml", "KitEmail", Format.GraphML),
        #GraphDescription("./input/LFR-1000.graph", "LFR-1000", Format.METIS),
        #GraphDescription("./input/PGPgiantcompo.graph", "PGP", Format.METIS),
        #GraphDescription("./input/karate.graph", "Karate", Format.METIS),
        #GraphDescription("./input/bter-graph-cc0.7-small.graph", "BTER", Format.EdgeListCommaOne),
        #GraphDescription("./input/ErdosRenyi.graph", "ErdosRenyi", Format.METIS),
        #GraphDescription("./input/jazz.graph", "Jazz", Format.METIS),

		GraphDescription("./input/cit-HepTh.edgelist-t0.graph", "HepTh", Format.EdgeList, continuous=False, separator='\t', firstNode=0),
		GraphDescription("./input/cit-HepPh.edgelist-t0.graph", "HepPh", Format.EdgeList, continuous=False, separator='\t', firstNode=0),
		GraphDescription("./input/soc-Epinions1.edgelist-t0.graph", "Epinions", Format.EdgeListTabZero),
		GraphDescription("./input/as20000102.txt", "AS", Format.EdgeList, continuous=False, separator='\t', firstNode=1)
        ]

#Returns a list of property calculation descriptions
def getProperties():
    return [
        P_General(),
        P_Ratios(),
        P_Community(),
        P_Diameter(),
        P_DegreeDistribution(),
        P_Centrality(),
        P_Components()
		#P_KolmogorowSmirnow()
        ]
