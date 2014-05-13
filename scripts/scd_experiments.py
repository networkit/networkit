from NetworKit import *
from selectivecommunity import *

import logging
import numpy
import os
import pandas
import math

def jaccard(A, B):
	i = len(A.intersection(B))
	u = len(A) + len(B) - i
	return i / u


def loadLFR(dir, n, mu, minc=50, maxc=250, k=15, maxk=100, t1=2.5, t2=1.5):
	mu = round(mu, 2)
	fpattern = "{0}_n={1}_mu={2}_minc={3}_maxc_{4}_k={5}_maxk={6}_t1={7}_t2={8}.dat"
	communityFileName = fpattern.format("community", n, mu, minc, maxc, k, maxk, t1, t2)
	logging.info("loading {0}".format(communityFileName))
	truth = community.readCommunities(os.path.join(dir, communityFileName), formats.edgelist_tab_zero)
	graphFileName = fpattern.format("network", n, mu, minc, maxc, k, maxk, t1, t2)
	logging.info("loading {0}".format(graphFileName))
	G = readGraph(os.path.join(dir, graphFileName), format=formats.edgelist_tab_one)
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


print("creating value ranges")
alphaRange = numpy.around(numpy.arange(0.01, 0.26, 0.01), 2)
epsilonRange = numpy.around(numpy.arange(1e-6, 10 * 1e-6, 1e-6), 6)

# available range of LFR parameters
muRange = numpy.around(numpy.arange(0.1, 1.0, 0.05), 2)

#nRange = [1000, 10000, 100000, 1000000]
nRange = [1000, 10000]

print("done value ranges")



def parameterStudyPageRankNibble(LFRDir, nSeeds=100):

	logging.info("starting parameter study")
	# range for alpha

	meanAccuracies = {}
	for n in nRange:
		for mu in muRange:
			print(mu)
			# load graph and ground truth
			(G, truth) = loadLFR(LFRDir, n, mu)
			seeds = [G.randomNode() for i in range(nSeeds)]	# reuse the same seeds

			# make data tables
			meanAccuracy = numpy.zeros((len(alphaRange), len(epsilonRange)))

			for alpha in alphaRange:
				for epsilon in epsilonRange:
					pageRankNibble = PageRankNibble(G, alpha, epsilon)
					(communities, accuracy) = LFRAccuracy(pageRankNibble, G, seeds, truth)
					# aggregate accuracy over seeds
					a, e = numpy.where(alphaRange == alpha), numpy.where(epsilonRange == epsilon)
					meanAccuracy[a,e] = numpy.mean([a for a in accuracy.values()])
					# print("alpha={0}, epsilon={1} -> maxAccuracy={3}, meanAccuracy={4}".format(alpha, epsilon, maxAccuracy, meanAccuracy))

					# make data frames
					# meanAccuracyDf = pandas.DataFrame(meanAccuracy, index=pandas.MultiIndex.from_product([alphaRange, epsilonRange], names=['alpha', 'epsilon']))
					# print(meanAccuracyDf)
			meanAccuracies[(n, mu)] = meanAccuracy

	return meanAccuracies
