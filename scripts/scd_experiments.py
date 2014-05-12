from NetworKit import *
from selectivecommunity import *

import logging
import numpy

def jaccard(A, B):
	i = len(A.intersection(B))
	u = len(A) + len(B) - i
	return i / u


def loadLFR(n, mu, minc=50, maxc=250, k=15, maxk=100, t1=2.5, t2=1.5):
	mu = round(mu, 2)
	fpattern = "{0}_n={1}_mu={2}_minc={3}_maxc_{4}_k={5}_maxk={6}_t1={7}_t2={8}.dat"
	communityFileName = fpattern.format("community", n, mu, minc, maxc, k, maxk, t1, t2)
	logging.info("loading {0}".format(communityFileName))
	truth = community.readCommunities(communityFileName, formats.edgelist_tab_zero)
	graphFileName = fpattern.format("network", n, mu, minc, maxc, k, maxk, t1, t2)
	logging.info("loading {0}".format(graphFileName))
	G = readGraph(graphFileName, format=formats.edgelist_tab_one)
	return (G, truth)



def LFRAccuracy(SCDAlgo, G, seeds, truth):

	SCDAlgo.run(seeds)
	communities = SCDAlgo.getResult()
	accuracy = {}

	for s in seeds:
		community = communities[s]
		trueCommunity = truth.getMembers(truth.subsetOf(s))
		accuracy[s] = jaccard(community, trueCommunity)

		logging.debug("seed {0}: jaccard {1}".format(s, accuracy[s]))
	return (communities, accuracy)





def parameterStudyPageRankNibble():
	logging.info("starting parameter study")
	# range for alpha
	alphaRange = numpy.arange(0.01, 0.26, 0.01)
	epsilonRange = numpy.arange(0.01, 0.11, 0.01)

	# available range of LFR parameters
	muRange = numpy.arange(0.1, 1.0, 0.05)
	#nRange = [1000, 10000, 100000, 1000000]
	nRange = [1000]

	# number of seeds
	nSeeds = 100

	for n in nRange:
		n = str(n)
		print(n)
		for mu in muRange:
			# load graph and ground truth
			(G, truth) = loadLFR(n, mu)
			seeds = [G.randomNode() for i in range(nSeeds)]	# reuse the same seeds

			for alpha in alphaRange:
				for epsilon in epsilonRange:
						alpha, epsilon = round(alpha, 2), round(epsilon, 2)
						pageRankNibble = PageRankNibble(G, alpha, epsilon)
						(communities, accuracy) = LFRAccuracy(pageRankNibble, G, seeds, truth)
						# aggregate accuracy over seeds
						minAccuracy, maxAccuracy = numpy.min([a for a in accuracy.values()]), numpy.max([a for a in accuracy.values()])
						meanAccuracy = numpy.mean([a for a in accuracy.values()])
						print("alpha={0}, epsilon={1} -> minAccuray={2}, maxAccuracy={3}, meanAccuracy={4}".format(alpha, epsilon, minAccuracy, maxAccuracy, meanAccuracy))

						
