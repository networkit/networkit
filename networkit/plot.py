from networkit import *
import matplotlib.pyplot as plt
import numpy as np
import operator

try:
	import seaborn
	seaborn.set_style("whitegrid")
except ImportError as importError:
	print("WARNING: module 'seaborn' is not installed, but recommended")
	print(importError)


def degreeDistribution(G):
	"""Plots the degree distribution of the given network."""
	dd = properties.degreeDistribution(G)
	plt.yscale("symlog")
	plt.xscale("log")
	plt.xlabel("nodes")
	plt.ylabel("degree")
	plt.plot(dd)

# TODO: hop plot

# TODO: core decomposition sequence

# TODO: clustering coefficient per degree

def clusteringPerDegree(G):
	degs = properties.degreeSequence(G)
	cc = properties.ClusteringCoefficient.exactLocal(G)
	pairs = sorted(zip(degs, cc), key=operator.itemgetter(0))
	x, y = zip(*pairs)
	plt.xlabel("degree")
	plt.ylabel("clustering coefficient")
	plt.scatter(x, y)
	seaborn.jointplot(np.array(x), np.array(y), kind="reg")
