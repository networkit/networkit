"""
This module handles compatibility between NetworKit and NetworkX
"""

import graph

# non standard library modules
try:
	import networkx as nx
except ImportError:
	print("""WARNING: module 'networkx' not installed, which is required by some
						functions.""")
########  CONVERSION ########

def nx2nk(nxG, weightAttr=None):
	""" 
	Convert a networkx.Graph to a NetworKit.Graph
		:param weightAttr: the edge attribute which should be treated as the edge weight
	 """
	# TODO: consider weights
	n = nxG.number_of_nodes()
	nkG = graph.Graph(n)
	
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
