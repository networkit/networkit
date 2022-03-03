"""
This module handles compatibility between NetworKit and NetworkX
"""

# local imports
from . import graph
import warnings
from .support import MissingDependencyError

# non standard library modules / external
try:
	import networkx as nx
except ImportError:
	have_nx = False
else:
	have_nx = True

########  CONVERSION ########

def nx2nk(nxG, weightAttr=None):
	"""
	nx2nk(nxG, weightAttr=None)

	Convert a networkx.Graph to a networkit.Graph.

	Parameters
	----------
	nxG : networkx.Graph
		The input networkx graph.
	weightAttr : str, optional
		The edge attribute which should be treated as the edge weight. Default: None
	"""

	if not have_nx:
		raise MissingDependencyError("networkx")
	# map networkx node ids to consecutive numerical node ids
	idmap = dict((id, u) for (id, u) in zip(nxG.nodes(), range(nxG.number_of_nodes())))
	z = max(idmap.values()) + 1
	# print("z = {0}".format(z))

	if weightAttr is not None:
		nkG = graph.Graph(z, weighted=True, directed=nxG.is_directed())
		for (u_, v_) in nxG.edges():
			u, v = idmap[u_], idmap[v_]
			w = nxG[u_][v_][weightAttr]
			nkG.addEdge(u, v, w)
	else:
		nkG = graph.Graph(z, directed=nxG.is_directed())
		for (u_, v_) in nxG.edges():
			u, v = idmap[u_], idmap[v_]
			assert (u < z)
			assert (v < z)
			nkG.addEdge(u, v)

	assert (nkG.numberOfNodes() == nxG.number_of_nodes())
	assert (nkG.numberOfEdges() == nxG.number_of_edges())
	return nkG

def nk2nx(nkG):
	""" 
	nk2nx(nkG)

	Convert a networKit.Graph to a networkx.Graph.

	Parameters
	----------
	nkG : networkit.Graph
		The input NetworKit graph.
	"""

	if not have_nx:
		raise MissingDependencyError("networkx")

	if nkG.isDirected():
		nxG = nx.DiGraph()
	else:
		nxG = nx.Graph()
	nxG.add_nodes_from(nkG.iterNodes())
	if nkG.isWeighted():
		for u, v, w in nkG.iterEdgesWeights():
			nxG.add_edge(u, v, weight=w)
	else:
		nxG.add_edges_from(nkG.iterEdges())

	assert (nkG.numberOfNodes() == nxG.number_of_nodes())
	assert (nkG.numberOfEdges() == nxG.number_of_edges())
	return nxG
