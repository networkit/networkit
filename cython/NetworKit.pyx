# cython: language_level=3

#includes


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

# standard library modules
import math
import textwrap

# non standard library modules
import networkx as nx
import tabulate


# local modules
import termgraph



# helper functions

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes


# Cython class definitions

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count) except +
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
		string toString()
		string getName()
		

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
		return self._this.hasEdge(u, v)
		
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

	def toString(self):
		return self._this.toString()

	def getName(self):
		return self._this.getName()
	
	
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
		

cdef extern from "../src/io/EdgeListIO.h":
	cdef cppclass _EdgeListIO "NetworKit::EdgeListIO":
		_EdgeListIO() except +
		_EdgeListIO(char separator, node firstNode) except +
		_Graph read(string path)
		void write(_Graph G, string path)

cdef class EdgeListIO:

	cdef _EdgeListIO _this

	def __cinit__(self, separator, firstNode):
		cdef char sep = stdstring(separator)[0]
		self._this = _EdgeListIO(sep, firstNode)

	def read(self, path):
		return Graph().setThis(self._this.read(stdstring(path)))

	def write(self, Graph G not None, path):
		self._this.write(G._this, stdstring(path))

# cdef extern from "../src/io/EdgeListIO.h":
# 	cdef cppclass _EdgeListIO "NetworKit::EdgeListIO":
# 		_EdgeListIO() except +
# 		_EdgeListIO(char separator, node firstNode) except +
# 		_Graph read(string path)
# 		void write(_Graph G, string path)

# cdef class EdgeListIO:

# 	cdef _EdgeListIO _this
	
# 	def __cinit__(self, string separator, node firstNode):
# 		cdef char sep = separator[0]
# 		self._this = _EdgeListIO(sep, firstNode)
		
# 	def read(self, path):
# 		return Graph().setThis(self._this.read(stdstring(path)))
	
# 	def write(self, Graph G not None, path):
# 		 # string needs to be converted to bytes, which are coerced to std::string
# 		self._this.write(G._this, stdstring(path))




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
		_Louvain(string par, double gamma)
		_Clustering run(_Graph _G)
		
cdef class Louvain(Clusterer):
	cdef _Louvain _this
	
	def __cinit__(self, par, gamma=1.0):
		self._this = _Louvain(stdstring(par), gamma)
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))



# this is an example for using static methods
cdef extern from "../src/properties/GraphProperties.h" namespace "NetworKit::GraphProperties":
	# static methods live in the class namespace, so declare them here
	pair[count, count] minMaxDegree(_Graph _G)
	double averageDegree(_Graph _G)
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
	def averageDegree(Graph G not None):
		return averageDegree(G._this)

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
		double getQuality(_Clustering _zeta, _Graph _G)
		

cdef class Modularity:
	cdef _Modularity _this
	
	def getQuality(self, Clustering zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)

# cdef extern from "../src/dcd/DynamicLabelPropagation.h":
# 	cdef cppclass _DynamicLabelPropagation "NetworKit::DynamicLabelPropagation":
# 		_DynamicLabelPropagation() except +
# 		_DynamicLabelPropagation(_Graph G, count theta, string strategyName) except +
# 		_Clustering run()
# 		string toString()
# 
# class DynamicCommunityDetector:
# 	pass	
# 		
# cdef class DynamicLabelPropagation:
# 	cdef _DynamicLabelPropagation _this
# 	
# 	def __cinit__(self, Graph G not None, theta, strategyName):
# 		self._this = _DynamicLabelPropagation(G._this, theta, stdstring(strategyName))
# 		
# 	def run(self):
# 		self._this.run()


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

def readGraph(path, format=None):
	"""	Read graph and return a NetworKit::Graph"""
	# TODO: detect file format by looking at the file content
	if format is None:	# look at file extension
		if path.endswith(".graph"):
			reader = METISGraphReader()
		else:
			raise Exception("unknown graph file format")
	else:
		if (format is "metis"):
			reader = METISGraphReader()
		elif (format is "edgelist"):
			reader = EdgeListIO('\t', 1)

	with open(path, "r") as file:	# catch a wrong path before it crashes the interpreter
		G = reader.read(path)
		return G

	return None
	



########  CONVERSION ########

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
	

########  PROPERTIES ########


def properties(nkG):
	print("converting to NetworX.Graph for some properties....")
	nxG = nxG = nk2nx(nkG)

	n = nkG.numberOfNodes()
	m = nkG.numberOfEdges()
	minMaxDeg = GraphProperties.minMaxDegree(nkG)
	avgDeg = GraphProperties.averageDegree(nkG)

	# calculating diameter for small graphs
	if (n < 100):
		print("calculating diameter...")
		dia = nx.diameter(nxG)
	else:
		dia = "[omitted]"


	# calculate eccentricity
	if (n < 100):
		eccentricities = nx.eccentricity(G)
		ecc = sum(val for val in eccentricities.values()) / n
	else:
		ecc = "[omitted]"



	# perform PLP and PLM community detection
	print("performing community detection: PLP")
	# TODO: avoid printout of status bar
	plp = LabelPropagation()
	zetaPLP = plp.run(nkG)
	ncomPLP = zetaPLP.numberOfClusters()
	modPLP = Modularity().getQuality(zetaPLP, nkG)
	print("performing community detection: PLM")
	PLM = Louvain("balanced")
	zetaPLM = PLM.run(nkG)
	ncomPLM = zetaPLM.numberOfClusters()
	modPLM = Modularity().getQuality(zetaPLM, nkG)

	# degree histogram
	
	print("calculating properties")	
	histo = nx.degree_histogram(nxG)
	(labels, histo) = compressHistogram(histo, nbins=25)

	# connected components
	components nx.connected_components(nxG)
	nComponents = len(components)
	componentSizes = [len(component) for component in components]
	componentSizes.sort(reverse=True) # sort in descending order
	sizeLargestComponent = componentsSizes[0]

	# betweenness centrality
	# TODO: average betweenness centrality?

	props = {
		 "name": nkG.getName(),
		 "n": n,
		 "m": m,
		 "minDeg": minMaxDeg[0],
		 "maxDeg": minMaxDeg[1],
		 "avgDeg": 
		 "avglcc": GraphProperties.averageLocalClusteringCoefficient(nkG),
		 "nComponents": nComponents,
		 "sizeLargestComponent": sizeLargestComponent,
		 "dia": dia,
		 "ecc": ecc,
		 "isolates": len(nx.isolates(nxG)),
		 "loops": len(nxG.selfloop_edges()),
		 "ncomPLP": ncomPLP,
		 "modPLP": modPLP,
		 "ncomPLM": ncomPLM,
		 "modPLM": modPLM,
		 "dens": nx.density(nxG),
		 "assort": nx.degree_assortativity_coefficient(nxG),
		 "cliques": len(list(nx.find_cliques(nxG))),
		 "histo": (labels, histo),
		 }

	return props


def showProperties(nkG):
	props = properties(nkG)
	basicProperties = [
		["nodes (n)", props["n"]],
		["edges (m)", props["m"]],
		["min. degree", props["minDeg"]],
		["max. degree", props["maxDeg"]],
		["avg. degree", props["avgDeg"]],
		["isolated nodes", props["isolates"]],
		["self-loops", props["loops"]],
		["density", "{0:.6f}".format(props["dens"])]
	]
	pathStructure = [
		["connected components", props["nComponents"]],
		["size of largest component", props["sizeLargestComponent"]],
		["diameter", props["dia"]],
		["avg. eccentricity", props["ecc"]],
	]
	
	miscProperties = [
		["degree assortativity", "{0:.6f}".format(props["assort"])],
		["cliques", props["cliques"]]
	]

	communityStructure = [
		["avg. local clustering coefficient", "", "{0:.6f}".format(props["avglcc"])],
		["PLP community detection", "", ""],
		["", "communities", props["ncomPLP"]],
		["", "modularity", "{0:.6f}".format(props["modPLP"])],
		["PLM community detection", "", ""],
		["", "communities", props["ncomPLM"]],
		["", "modularity", "{0:.6f}".format(props["modPLM"])],

	]

	print()
	print("Basic Properties")
	print(tabulate.tabulate(basicProperties))
	print("Path Structure")
	print(tabulate.tabulate(pathStructure))
	print("Miscellaneous")
	print(tabulate.tabulate(miscProperties))
	print("Community Structure")
	print(tabulate.tabulate(communityStructure))
	print("Degree Distribution")
	print("-------------------")
	(labels, histo) = props["histo"]
	termgraph.graph(labels, histo)



def showPropertiesOld(nkG):

	propertiesTextBlock = """
	Graph Properties: {name}
	========================

	Basic Properties
	----------------
	- nodes (n){n:16}
	- edges (m)			{m}
	- min. degree 			{minDeg}
	- max. degree 			{maxDeg}
	- isolated nodes 		{isolates}
	- self-loops			{loops}
	- density			{dens:10.6f}


	Path Structure
	--------------
	- connected components 		{components}
	- diameter			{dia}

	Miscellaneous
	-------------
	- degree assortativity			{assort:10.6f}
	- cliques				{cliques}

	Community Structure
	-------------------
	- avg. local clustering coefficient 		{avglcc:10.6f}
	- PLP community detection
		- communities				{ncomPLP}
		- modularity 			{modPLP:10.6f}
	- PLM community detection
		- communities			{ncomPLM}
		- modularity 			{modPLM:10.6f}


	Distributions
	-------------

	- degree distribution

	"""

	props = properties(nkG)
	print(textwrap.dedent(propertiesTextBlock.format(**props)))
	(labels, histo) = props["histo"]
	termgraph.graph(labels, histo)

	

def compressHistogram(hist, nbins=20):
	""" Compress a histogram to a number of bins"""
	compressed = [None for i in range(nbins)]
	labels = [None for i in range(nbins)]
	nbinsprev = len(hist)
	binwidth = math.ceil(nbinsprev / nbins)
	for i in range(nbins):
		l = i * binwidth
		u = (i+1) * binwidth
		compressed[i] = sum(hist[l : u])
		labels[i] = "{0}-".format(l)
	return (labels, compressed)
		

	
def printDegreeHistogram(nxG, nbins=25):
	""" Prints a degree histogram as a bar chart to the terminal"""
	hist = nx.degree_histogram(nxG)
	
	(labels, hist) = compressHistogram(hist, nbins)
	termgraph.graph(labels, hist)
	


def hpProperties(nkG):
	""" For large graphs: get an overview of some properties"""
	print("min/max degree:")
	
	
	
# NetworKit algorithm engineering workflows

# class DynamicCommunityDetectionWorkflow:
# 	
# 	def __init__(self):
# 		clusterings = []	# list of clusterings
# 	
# 	
# 	def start(self, nMax, deltaT):
# 		
# 		self.G = Graph(0)
# 		self.Gproxy = GraphEventProxy(self.G)
# 		#self.generator = DynamicBarabasiAlbertGenerator(self.Gproxy)
# 		self.dcd = DynamicLabelPropagation(self.Gproxy)
# 		
# 		while (self.G.numberOfNodes() < nMax):
# 			self.generator.generate()
# 			if (self.G.time() % deltaT) == 0:
# 				zeta = self.dcd.run()
# 				self.clusterings.append(zeta)
	
class GraphConverter:
	
	def __init__(self, reader, writer):
		self.reader = reader
		self.writer = writer
		
	def convert(self, inPath, outPath):
		G = self.reader.read(inPath)
		self.writer.write(G, outPath)
		
	def __str__(self):
		return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
	
	readers = {"metis": METISGraphReader, "edgelist" : EdgeListIO}	
	writers = {"edgelist": EdgeListIO}
	
	reader = readers[fromFormat]()
	writer = writers[toFormat]()
	
	return GraphConverter(reader, writer)

	