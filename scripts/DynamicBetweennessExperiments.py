from NetworKit import *
from dynamic import *


def removeAndAddEdges(G, nEdges, batchSize=1, tabu=None):
	if nEdges > G.numberOfEdges():
		raise Error("G does not have enough edges")

	# select random edges for removal
	removed = set()
	while len(removed) < nEdges:
		(u, v) = G.randomEdge()
		if not tabu.hasEdge(u, v):	# exclude all edges in the tabu graph
			removed.add((u, v))

	# build event streams
	removeStream = []
	for (u, v) in removed:
		removeStream.append(GraphEvent(GraphEvent.EDGE_REMOVAL, u, v, 0))
	addStream = []
	i = 0
	for (u, v) in removed:
		if (i % batchSize) == 0:
			addStream.append(GraphEvent(GraphEvent.TIME_STEP, 0, 0, 0))
		addStream.append(GraphEvent(GraphEvent.EDGE_ADDITION, u, v, 1.0))
		i += 1

	return (removeStream, addStream)



def test():
	G = readGraph("../input/PGPgiantcompo.graph")
	T = graph.SpanningForest(G).generate()

	(removeStream, addStream) = removeAndAddEdges(G, 100, batchSize=5, tabu=T)
	return (removeStream, addStream)
