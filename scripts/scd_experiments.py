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
	print("loading {0}".format(communityFileName))
	truth = community.readCommunities(os.path.join(dir, communityFileName), formats.edgelist_tab_one)
	graphFileName = fpattern.format("network", n, mu, minc, maxc, k, maxk, t1, t2)
	print("loading {0}".format(graphFileName))
	G = readGraph(os.path.join(dir, graphFileName), format=formats.edgelist_tab_one)
	return (G, truth)



def LFRAccuracy(SCDAlgo, G, seeds, truth):
	timer = stopwatch.Timer()
	#
	SCDAlgo.run(seeds)
	#
	time = timer.stop()
	communities = SCDAlgo.getResult()
	accuracy = {}

	for s in seeds:
		community = communities[s]
		trueCommunity = truth.getMembers(truth.subsetOf(s))
		accuracy[s] = jaccard(community, trueCommunity)

		logging.debug("seed {0}: jaccard {1}".format(s, accuracy[s]))
	return (communities, accuracy, time)


class PseudoSCD:

	def __init__(self, G, cdAlgo):
		self.cdAlgo = cdAlgo
		self.zeta = None
		self.communities = {}
		self.G = G

	def run(self, seeds):
		if not self.zeta:
			self.zeta = self.cdAlgo.run(self.G)
		for s in seeds:
			self.communities[s] = self.zeta.getMembers(self.zeta.subsetOf(s))

	def getResult(self):
		return self.communities


def testAccuracyOnLFR(AlgoClass, algoParams, LFRDir, nSeeds=100):

	print("start LFR test run")
	# range for alpha

	meanAccuracy = {}
	timings = []
	for n in nRange:
		for mu in muRange:
			print(mu)
			# load graph and ground truth
			(G, truth) = loadLFR(LFRDir, n, mu)
			seeds = [G.randomNode() for i in range(nSeeds)]	# reuse the same seeds
			algo = AlgoClass(G, *algoParams)
			(communities, accuracy, time) = LFRAccuracy(algo, G, seeds, truth)
			timings.append(time)

			#print("accuracy values: ", [a for a in accuracy.values()])
			# aggregate accuracy over seeds
			# store per LFR params
			meanAccuracy[(n, mu)] = numpy.mean([a for a in accuracy.values()])
	print("done")
	return (meanAccuracy, timings)




# available range of LFR parameters
muRange = numpy.around(numpy.arange(0.0, 1.0, 0.05), 2)

#nRange = [1000, 10000, 100000, 1000000]
nRange = [10000]

# parameter ranges for PageRankNibble

alphaRange = numpy.around(numpy.arange(0.01, 0.26, 0.01), 2)
epsilonRange = numpy.around(numpy.arange(1e-5, 10 * 1e-5, 1e-5), 5)

def parameterStudyPageRankNibble(LFRDir, nSeeds=100):

	print("starting parameter study")
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

			runs = 0
			print("running experiments")
			for alpha in alphaRange:
				for epsilon in epsilonRange:
					pageRankNibble = PageRankNibble(G, alpha, epsilon)
					(communities, accuracy, time) = LFRAccuracy(pageRankNibble, G, seeds, truth)
					runs += 1
					#print("accuracy values: ", [a for a in accuracy.values()])
					# aggregate accuracy over seeds
					a, e = numpy.where(alphaRange == alpha), numpy.where(epsilonRange == epsilon)
					meanAccuracy[a,e] = numpy.mean([a for a in accuracy.values()])
					# print("alpha={0}, epsilon={1} -> maxAccuracy={3}, meanAccuracy={4}".format(alpha, epsilon, maxAccuracy, meanAccuracy))

					# make data frames
					# meanAccuracyDf = pandas.DataFrame(meanAccuracy, index=pandas.MultiIndex.from_product([alphaRange, epsilonRange], names=['alpha', 'epsilon']))
					# print(meanAccuracyDf)
			meanAccuracies[(n, mu)] = meanAccuracy

	print("done")
	return meanAccuracies


# parameter ranges for SelSCAN

kappaRange = range(1, 4)
epsilonRange = numpy.arange(0.0, 1.0, 0.1)

def parameterStudySelSCAN(LFRDir, nSeeds=100):

	print("starting parameter study")
	# range for alpha

	meanAccuracies = {}
	for n in nRange:
		for mu in muRange:
			print(mu)
			# load graph and ground truth
			(G, truth) = loadLFR(LFRDir, n, mu)
			seeds = [G.randomNode() for i in range(nSeeds)]	# reuse the same seeds

			# make data tables
			meanAccuracy = numpy.zeros((len(epsilonRange), len(kappaRange)))

			runs = 0
			print("running experiments")
			for epsilon in epsilonRange:
				for kappa in kappaRange:
					selSCAN = SelSCAN(G, kappa, epsilon)
					(communities, accuracy, time) = LFRAccuracy(selSCAN, G, seeds, truth)
					runs += 1
					#print("accuracy values: ", [a for a in accuracy.values()])
					# aggregate accuracy over seeds
					a, e = numpy.where(epsilonRange == epsilon), numpy.where(kappaRange == kappa)
					meanAccuracy[a,e] = numpy.mean([a for a in accuracy.values()])

			meanAccuracies[(n, mu)] = meanAccuracy

	print("done")
	return meanAccuracies
