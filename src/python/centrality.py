from _NetworKit import Betweenness, Betweenness2, PageRank, EigenvectorCentrality, DegreeCentrality, ApproxBetweenness, ApproxBetweenness2


def ranking(G, algorithm=Betweenness, normalized=False):
	""" Return a ranking of nodes by the specified centrality type"""
	# FIXME: some centrality algorithms take more parameters
	centrality = algorithm(G, normalized)
	centrality.run()
	return centrality.ranking()

def scores(G, algorithm=Betweenness, normalized=False):
	""" Return the centrality scores of nodes using the specified centrality type"""
	centrality = algorithm(G, normalized)
	centrality.run()
	return centrality.scores()
