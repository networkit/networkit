from NetworKit import *


def jaccard(A, B):
	i = len(A.intersection(B))
	u = len(A) + len(B) - i
	return i / u


def loadLFR(n, mu, minc, maxc, k, maxk, t1, t2):
	fname = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}"
	truth = community.readCommunities(fname.format("community", mu, minc, maxc, k, maxk, t1, t2), formats.edgelist_tab_zero)
	G = readGraph(fname.format("network", mu, minc, maxc, k, maxk, t1, t2))
	return (G, truth)



def LFRBenchmark(SCDAlgo, G, truth):

	G = readGraph(graphPath)
	communities = SCDAlgo.run(seeds)


	for s in seed:
		community = communities[s]
		trueCommunity = truth.getMembers(truth.getSubset(s))

		j = jaccard(community, trueCommunity)

	return j
