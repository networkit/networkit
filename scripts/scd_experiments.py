from NetworKit import *
from selectivecommunity import *

import numpy

def jaccard(A, B):
	i = len(A.intersection(B))
	u = len(A) + len(B) - i
	return i / u


def loadLFR(n, mu, minc=50, maxc=250, k=15, maxk=100, t1=2.5, t2=1.5):
	fpattern = "{0}_n={1}_mu={2}_minc={3}_maxc={4}_k={5}_maxk={6}_t1={7}_t2={8}.dat"
	fname = fpattern.format("community", n, mu, minc, maxc, k, maxk, t1, t2)
	print("loading ", fname)
	truth = community.readCommunities(fname, formats.edgelist_tab_zero)
	G = readGraph(fname.format("network", mu, minc, maxc, k, maxk, t1, t2))
	return (G, truth)



def LFRAccuracy(SCDAlgo, G, seeds, truth):

	SCDAlgo.run(seeds)
	communities = SCDAlgo.getResult()
	accuracy = {}

	for s in seeds:
		community = communities[s]
		trueCommunity = truth.getMembers(truth.subsetOf(s))
		accuracy[s] = jaccard(community, trueCommunity)

		print("seed {0}: jaccard {1}".format(s, j))
	return (communities, accuracy)





def parameterStudyPageRankNibble():
	print("starting parameter study")
	# range for alpha
	alphaRange = numpy.arange(0.01, 0.26, 0.01)
	epsilonRange = numpy.arange(0.01, 0.11, 0.01)

	# available range of LFR parameters
	muRange = numpy.arange(0.1, 1.0, 0.05)
	nRange = [1000, 10000, 100000, 1000000]

	# number of seeds
	nSeeds = 10

	for n in nRange:
		n = str(n)
		print(n)
		for mu in muRange:
			print(n, mu)
			# load graph and ground truth
			(G, truth) = loadLFR(n, mu)

			for alpha in alphaRange:
				for epsilon in epsilonRange:
						pageRankNibble = PageRankNibble(G, alpha, epsilon)
						seeds = [G.randomNode() for i in range(nSeeds)]
						accuracy = LFRAccuracy(pageRankNibble, G, seeds, truth)
						# aggregate accuracy over seeds
						minAccuracy, maxAccuracy = numpy.min(accuracy.values()), maxAccuracy
						meanAccuracy = numpy.mean(accuracy.values())
						print("alpha={0}, epsilon={1} -> {2}, {3}, {4}".format(alpha, epsilon, minAccuracy, maxAccuracy, meanAccuracy))
