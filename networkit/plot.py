from networkit import *
import matplotlib.pyplot as plt
import numpy as np
import operator
import pandas

try:
	import seaborn
	seaborn.set_style("whitegrid")
except ImportError as importError:
	print("WARNING: module 'seaborn' is not installed, plotting functionality will be limited")
	print(importError)


def degreeDistribution(G, **kwargs):
	"""Plots the degree distribution of the given network."""
	dd = properties.degreeDistribution(G)
	plt.yscale("symlog")
	plt.xscale("log")
	plt.xlabel("nodes")
	plt.ylabel("degree")
	plt.plot(dd)

# TODO: hop plot

def coreDecompositionSequence(G, **kwargs):
	""" Plots the core decomposition sequence of G, i.e. the size of the k-shell for the core number k"""
	shells = properties.CoreDecomposition(G).run().shells()
	data = pandas.DataFrame({"k": range(len(shells)), "n_k": [len(shell) for shell in shells]})
	plt.xlabel("core number")
	plt.ylabel("number of nodes")
	plt.plot(data["k"], data["n_k"], **kwargs)


def clusteringPerDegree(G, **kwargs):
	""" Plots the local clustering coefficient for nodes with specific degree"""
	degs = properties.degreeSequence(G)
	cc = properties.ClusteringCoefficient.exactLocal(G)
	data = pandas.DataFrame({"deg": degs, "cc" : cc})
	data = data.groupby("deg", as_index=False).mean()
	jointplot = seaborn.jointplot("deg", "cc", data, kind="reg", ylim=(0, 1), **kwargs)
