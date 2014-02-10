from _NetworKit import Graph, GraphEvent, DynamicCommunityDetection, DGSStreamParser, GraphUpdater


def graphFromStream(stream):
	G = Graph()
	gu = GraphUpdater(G)
	gu.update(stream)
	return G