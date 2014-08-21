'''
Created on 09.04.2013

@author: cls
'''

# external imports
import random
import networkx
import itertools
import math


class GraphGenerator:
	
	def makeWeightedErdosRenyiGraph(self, n=100, p=0.2, weightFunction=random.random):
		"""
		Make a weighted Erdos-Renyi random graph.

		:rtype G: The graph object
		"""
		G = networkx.generators.erdos_renyi_graph(n, p)
		for (u, v) in G.edges():
			G[u][v]["weight"] = weightFunction()
		return G
		
	def makeWeightedCircularGraph(self, n=100, p=0.2, weightFunction=random.random):
		""" 
		Makes a circular graph with randomly weighted edges

		:rtype G: The graph object
		"""
		G = networkx.Graph()
		for u in range(n):
			G.add_node(u)
		for u in range(n):
			v = (u + 1) % n		
			G.add_edge(u, v)
			G[u][v]["weight"] = weightFunction()
			
		return G


	def makeClusteredRandomGraph(self, n, k, pin=0.5, pout=0.1):
		""" Makes a clustered random graph with intra-cluster edge probability pout + pin,
		inter-cluster edge probability pout,
		number of clusters k.

		:rtype G: The graph object
		"""
		
		def chunks(l, size):
			""" Yield successive n-sized chunks from l. """
			for i in range(0, len(l), size):
				yield l[i:i+size]
			
		
		G = networkx.generators.erdos_renyi_graph(n, pout)
		
		size = math.ceil(n / k)
		for c in chunks(G.nodes(), size):
			for (u,v) in itertools.combinations(c, 2):
				if (random.random() <= pin):
					G.add_edge(u, v)
		return G
		
