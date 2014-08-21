# external imports
import networkx as nx
import itertools


class GraphConstruct:
	"""
	This class is for the convenient construction of graphs 
	by chained method calls.
	"""
	
	def __init__(self):
		self.G = nx.Graph()
		self.vMax = 0
		
	def makeNode(self):
		v = self.vMax
		self.G.add_node(v)
		self.vMax += 1
		return v
	
	def addNode(self):
		self.makeNode()
		return self
	
	def addEdge(self, u, v):
		self.G.add_edge(u, v)
		return self
		
	def addClique(self, size, nodes=None):
		if nodes is None:
			nodes = [self.makeNode() for i in range(size)]
		for (u, v) in itertools.combinations(nodes, 2):
			self.G.add_edge(u, v)
		return self
			
	def getGraph(self):
		return self.G


if __name__ == "__main__":
	# example
	drawer.draw(GraphConstruct()
			.addClique(4)
			.addClique(4)
			.addClique(4)
			.addClique(4)
			.addEdge(0,4)
			.addEdge(8,12)
			.addEdge(9,15)
			.addEdge(3,9)
			.getGraph())
	savefig("/Users/cls/Downloads/conductance-example_2.pdf")
