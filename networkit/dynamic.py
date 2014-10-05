# extension imports
from _NetworKit import Graph, GraphEvent, DGSStreamParser, GraphUpdater


def graphFromStream(stream):
	""" Convenience function for creating a new graph from a stream of graph events

	Parameters
	----------
	stream : list of GraphEvent
		event stream
	"""
	G = Graph()
	gu = GraphUpdater(G)
	gu.update(stream)
	return G
