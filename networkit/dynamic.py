# extension imports
from .graph import Graph
from .dynamics import graphFromStream as dynamicsGraphFromStream

def graphFromStream(stream, weighted, directed):
	""" 
	graphFromStream(stream, weighted, directed)
	
	DEPRECATED. Use networkit.dynamics.graphFromStream instead.

	Convenience function for creating a new graph from a stream of graph events

	Parameters
	----------
	stream : list(networkit.dynamics.GraphEvent)
		Event stream
	weighted : bool
		Produce a weighted or unweighted graph
	directed : bool
		Produce a directed or undirected graph
	"""
	from warnings import warn
	warn("networkit.dynamic.graphFromStream is deprecated, use networkit.dynamics.graphFromStream")
	G = dynamicsGraphFromStream(stream, weighted, directed)
	return G
