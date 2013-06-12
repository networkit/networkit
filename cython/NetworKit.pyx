# type imports
from libc.stdint cimport uint64_t

# the C++ standard library
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

# NetworKit typedefs
ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef index cluster
ctypedef double edgeweight

# python module imports
import networkx as nx
# TODO: how to import termgraph



# helper functions

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes


# Cython class definitions

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(int) except +
		count numberOfNodes()
		count numberOfEdges()
		node addNode()
		void removeNode(node u)
		void addEdge(node u, node v, edgeweight w)
		# TODO: optional weight argument
		void removeEdge(node u, node v)
		bool hasEdge(node u, node v)
		edgeweight weight(node u, node v)
		vector[node] nodes()
		vector[pair[node, node]] edges()
		void markAsWeighted()
		bool isMarkedAsWeighted()
		

cdef class Graph:
	cdef _Graph _this
	
	def __cinit__(self, n=None):
		if n is not None:
			self._this = _Graph(n)
		
	# any _thisect which appears as a return type needs to implement setThis
	cdef setThis(self, _Graph other):
		#del self._this
		self._this = other
		return self
	
	def numberOfNodes(self):
		return self._this.numberOfNodes()
	
	def numberOfEdges(self):
		return self._this.numberOfEdges()
	
	def addNode(self):
		return self._this.addNode()
	
	def removeNode(self, u):
		self._this.removeNode(u)
		
	def addEdge(self, u, v, w=1.0):
		self._this.addEdge(u, v, w)
		
	def removeEdge(self, u, v):
		self._this.removeEdge(u, v)
		
	def hasEdge(self, u, v):
		self._this.hasEdge(u, v)
		
	def weight(self, u, v):
		return self._this.weight(u, v)
		
	def nodes(self):
		return self._this.nodes()
	
	def edges(self):
		return self._this.edges()
	
	def markAsWeighted(self):
		self._this.markAsWeighted()
	
	def isMarkedAsWeighted(self):
		return self._this.isMarkedAsWeighted()
	
	
cdef extern from "../src/graph/GraphGenerator.h":
	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
		_GraphGenerator() except +
		_Graph makeRandomGraph(count n, double p)


cdef class GraphGenerator:
	cdef _GraphGenerator _this
	
	def __cinit__(self):
		self._this = _GraphGenerator()
		
	
	def makeRandomGraph(self, n, p):
		cdef _Graph _G = self._this.makeRandomGraph(n, p)
		return Graph(0).setThis(_G)



cdef extern from "../src/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader":
		_METISGraphReader() except +
		_Graph read(string path)

cdef class METISGraphReader:
	cdef _METISGraphReader _this
	
	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		cdef _Graph _G = self._this.read(pathbytes)
		return Graph(0).setThis(_G)
	

cdef extern from "../src/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path)


cdef class DotGraphWriter:
	cdef _DotGraphWriter _this
	
	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		self._this.write(G._this, stdstring(path))
	

cdef extern from "../src/clustering/Clustering.h":
	cdef cppclass _Clustering "NetworKit::Clustering":
		_Clustering() except +
		count numberOfClusters()

cdef class Clustering:
	cdef _Clustering _this
	
	cdef setThis(self, _Clustering other):
		self._this = other
		return self

	def numberOfClusters(self):
		return self._this.numberOfClusters()

cdef class Clusterer:
	""" Abstract base class for static community detection algorithms"""
	def run(self, Graph G not None):
		raise NotImplementedError("abstract method")

cdef extern from "../src/community/LabelPropagation.h":
	cdef cppclass _LabelPropagation "NetworKit::LabelPropagation":
		_LabelPropagation() except +
		_Clustering run(_Graph _G)


cdef class LabelPropagation(Clusterer):
	cdef _LabelPropagation _this
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))
	
	

cdef extern from "../src/community/Louvain.h":
	cdef cppclass _Louvain "NetworKit::Louvain":
		_Louvain() except +
		_Clustering run(_Graph _G)
		
cdef class Louvain(Clusterer):
	cdef _Louvain _this
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))



# this is an example for using static methods
cdef extern from "../src/properties/GraphProperties.h" namespace "NetworKit::GraphProperties":
	# static methods live in the class namespace, so declare them here
	pair[count, count] minMaxDegree(_Graph _G)
	vector[count] degreeDistribution(_Graph _G)
	vector[double] localClusteringCoefficients(_Graph _G)
	double averageLocalClusteringCoefficient(_Graph _G)
	vector[double] localClusteringCoefficientPerDegree(_Graph _G)
	
	cdef cppclass _GraphProperties "NetworKit::GraphProperties":
		pass

cdef class GraphProperties:

	@staticmethod
	def minMaxDegree(Graph G not None):
		return minMaxDegree(G._this)

	@staticmethod
	def degreeDistribution(Graph G not None):
		return degreeDistribution(G._this)

	@staticmethod
	def averageLocalClusteringCoefficient(Graph G not None):
		return averageLocalClusteringCoefficient(G._this)


cdef extern from "../src/dynamics/GraphEventHandler.h":
	cdef cppclass _GraphEventHandler "NetworKit::GraphEventHandler":
		_GraphEventHandler() except +
		void onNodeAddition(node u)
		void onNodeRemoval(node u)
		void onEdgeAddition(node u, node v)
		void onEdgeRemoval(node u, node v)
		void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew)
		void onTimeStep()


cdef extern from "../src/dynamics/GraphEventProxy.h":
	cdef cppclass _GraphEventProxy "NetworKit::GraphEventProxy":
		_GraphEventProxy()	# nullary constructor not valid
		_GraphEventProxy(_Graph _G)
		void registerObserver(_GraphEventHandler* _observer)
		node addNode()
		void removeNode(node u)
		void addEdge(node u, node v)
		# TODO: optional edge weight
		void removeEdge(node u, node v)
		void setWeight(node u, node v, edgeweight w)
		void timeStep()
		
cdef class GraphEventProxy:
	cdef _GraphEventProxy _this
	
	def __cinit__(self, Graph G not None):
		self._this = _GraphEventProxy(G._this)
	# TODO: delegates



cdef extern from "../src/io/DGSReader.h":
	cdef cppclass _DGSReader "NetworKit::DGSReader":
		_DGSReader() except +
		void read(string path, _GraphEventProxy _proxy)
		
		
cdef class DGSReader:
	cdef _DGSReader _this
	
	def read(self, path, GraphEventProxy proxy not None):
		self._this.read(stdstring(path), proxy._this)
		
cdef extern from "../src/clustering/Modularity.h":
	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(const _Clustering _zeta, const _Graph _G)
		
		
cdef class Modularity:
	cdef _Modularity _this
	
	def getQuality(self, Clustering zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)

cdef extern from "../src/community/DynamicLabelPropagation.h":
	cdef cppclass _DynamicLabelPropagation "NetworKit::DynamicLabelPropagation":
		_DynamicLabelPropagation() except +
		_DynamicLabelPropagation(_Graph G, count theta, string strategyName) except +
		_Clustering run()
		string toString()

class DynamicCommunityDetector:
	pass	
		
cdef class DynamicLabelPropagation:
	cdef _DynamicLabelPropagation _this
	
	def __cinit__(self, Graph G not None, theta, strategyName):
		self._this = _DynamicLabelPropagation(G._this, theta, stdstring(strategyName))
		
	def run(self):
		self._this.run()


# FIXME:
# cdef extern from "../src/generators/DynamicBarabasiAlbertGenerator.h":
# 	cdef cppclass _DynamicBarabasiAlbertGenerator "NetworKit::DynamicBarabasiAlbertGenerator":
# 		_DynamicBarabasiAlbertGenerator() except +
# 		_DynamicBarabasiAlbertGenerator(_GraphEventProxy _Gproxy, count k) except +
# 		void initializeGraph()
# 		void generate()
# 		
# cdef class DynamicBarabasiAlbertGenerator:
# 	cdef _DynamicBarabasiAlbertGenerator _this
# 	
# 	def __cinit__(self, GraphEventProxy Gproxy not None, k):
# 		self._this = _DynamicBarabasiAlbertGenerator(Gproxy._this, k)
# 		
# 	def generate(self):
# 		self._this.generate()
		


# under construction

cdef extern from "../src/generators/PubWebGenerator.h":
	cdef cppclass _PubWebGenerator "NetworKit::PubWebGenerator":
		_PubWebGenerator() except +
		_PubWebGenerator(count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors) except +
		_Graph generate()
		
cdef class PubWebGenerator:
	cdef _PubWebGenerator _this
	
	def __cinit__(self, count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors):
		self._this = _PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)
		
	def generate(self):
		return Graph().setThis(self._this.generate())	

cdef extern from "../src/viz/ForceDirected.h":
	cdef cppclass _ForceDirected "NetworKit::ForceDirected":
		_ForceDirected() except +
		void draw(_Graph _G)
	
cdef class ForceDirected:
	cdef _ForceDirected _this
	
	def draw(self, Graph G not None):
		pass
	


# TODO: initialize log4cxx

#--------- NetworKit Python Shell functions ----------------#

def readGraph(path):
	"""	Automatically detect input format and read graph"""
	# TODO: detect file format by looking at the file content
	if path.endswith(".graph"):
		reader = METISGraphReader()
	else:
		raise Exception("unknown graph file format")
	G = reader.read(path)
	return G




def nx2nk(nxG, weightAttr=None):
	""" 
	Convert a networkx.Graph to a NetworKit.Graph
		:param weightAttr: the edge attribute which should be treated as the edge weight
	 """
	# TODO: consider weights
	n = nxG.number_of_nodes()
	cdef Graph nkG = Graph(n)
	
	if weightAttr is not None:
		nkG.markAsWeighted()
		for (u, v) in nxG.edges():
			w = nxG.edge[u][v][weightAttr]
			nkG.addEdge(u, v, w)
	else:
		for (u, v) in nxG.edges():
			nkG.addEdge(u, v)
	
	return nkG


def nk2nx(nkG):
	""" Convert a NetworKit.Graph to a networkx.Graph """
	nxG = nx.Graph()
	if nkG.isMarkedAsWeighted():
		for (u, v) in nkG.edges():
			nxG.add_edge(u, v, weight=nkG.weight(u, v))
	else:
		for (u, v) in nkG.edges():
			nxG.add_edge(u, v)
	return nxG
	
	
def properties(nkG):
	""" Get an overview of the properties for the graph"""
	nxG = nk2nx(nkG)
	if (nxG.markAsWeighted()):
		print("The graph marked as weighted")
	else: 
		print("The graph has not been marked as weighed")
	(n, m) = (nxG.number_of_nodes(), nxG.number_of_edges())
	print("Number of nodes \t n \t {0}".format(n))
	print("Number of edges \t n \t {0}".format(m))
	
	print("Degree distribution:")
	printDegreeHistogram(nxG)

	minMax = GraphProperties.minMaxDegree(nkG);
	print("Minimum degree \t {0}".format(minMax[0]))
	print("Maximum degree \t {0}".format(minMax[1]))

	degree_distribution = GraphProperties.degreeDistribution(nxG)
	print("Degree distribution:\n{0}".format(degree_distribution))
	
	isolated_nodes = nx.isolate(nxG)
	print("Isolated nodes: \t{0}".format(len(isolated_nodes))) 

	self_loops = nxG.selfloop_edges()
	print("Self-loop edges: \t{0}".format(len(self_loops)))

	print("Connected components: \t{0}".format(nx.number_connected_components(nxG)))

	average_Local_clustering = GraphProperties.averageLocalClusteringCoefficient
	print("Average local clustering coefficient: \t{0}".format(average_Local_clustering))

	plp = LabelPropagation()
	zeta = LabelPropagation.run(plp)
	n_clusters = zeta.numberOfClusters();
	print("Label propagation found: \t{0}".format(n_clusters))
	
	
	
def printDegreeHistogram(nxG):
	""" Prints a degree histogram as a bar chart to the terminal"""
	hist = nx.degree_histogram(nxG)
	labels = range(len(hist))
	# FIXME: termgraph.graph(labels, hist)
	


def hpProperties(nkG):
	""" For large graphs: get an overview of some properties"""
	print("min/max degree:")
	
	
def compactDegreeHistogram(nxG, nbins=10):
	"""
	Create a compact histogram for the degree distribution. 
	"""
	hist = nx.degree_histogram(nxG)
	maxDeg = len(hist) - 1
	binsize = len(hist) / nbins
	labels = []
	values = []
	
	# TODO: range(0, len(hist), binsize)
	
	return (labels, values)
	
	
# NetworKit algorithm engineering workflows

class DynamicCommunityDetectionWorkflow:
	
	def __init__(self):
		clusterings = []	# list of clusterings
	
	
	def start(self, nMax, deltaT):
		
		self.G = Graph(0)
		self.Gproxy = GraphEventProxy(self.G)
		#self.generator = DynamicBarabasiAlbertGenerator(self.Gproxy)
		self.dcd = DynamicLabelPropagation(self.Gproxy)
		
		while (self.G.numberOfNodes() < nMax):
			self.generator.generate()
			if (self.G.time() % deltaT) == 0:
				zeta = self.dcd.run()
				self.clusterings.append(zeta)
		
			