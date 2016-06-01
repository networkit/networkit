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


def nodeProperty(data, label, sorted=True, yscale="linear", xscale="linear"):
	""" General plotting function for a node property"""
	plt.yscale(yscale)
	plt.xscale(xscale)
	plt.xlabel("nodes")
	plt.ylabel(label)
	if sorted:
		data.sort(reverse=True)
	plt.plot(data)


def degreeDistribution(G, **kwargs):
	"""Plots the degree distribution of the given network."""
	dd = properties.degreeDistribution(G)
	plt.yscale("symlog")
	plt.xscale("log")
	plt.xlabel("nodes")
	plt.ylabel("degree")
	plt.plot(dd)

def connectedComponentsSizes(G, **kwargs):
	""" Plot the size distribution of connected components as a pie chart """
	csizes = components.ConnectedComponents(G).run().getComponentSizes()
	colors = seaborn.color_palette("Set2", 10)
	data = list(csizes.values())
	# explode the largest component pie piece
	maxi = data.index(max(data))
	explode = [0 for i in range(len(data))]
	explode[maxi] = 0.1
	# plot
	plt.figure(figsize=(5,5))
	plt.pie(data, colors=colors, autopct='%1.1f%%', explode=explode)

# TODO: hop plot

def coreDecompositionSequence(G, **kwargs):
	""" Plots the core decomposition sequence of G, i.e. the size of the k-shell for the core number k"""
	shells = centrality.CoreDecomposition(G).run().shells()
	data = pandas.DataFrame({"k": range(len(shells)), "n_k": [len(shell) for shell in shells]})
	plt.xlabel("core number")
	plt.ylabel("number of nodes")
	plt.plot(data["k"], data["n_k"], **kwargs)


def clusteringPerDegree(G, **kwargs):
	""" Plots the local clustering coefficient for nodes with specific degree"""
	degs = centrality.DegreeCentrality(G).run().scores()
	cc = centrality.LocalClusteringCoefficient(G).run().scores()
	data = pandas.DataFrame({"deg": degs, "cc" : cc})
	data = data.groupby("deg", as_index=False).mean()
	jointplot = seaborn.jointplot("deg", "cc", data, kind="reg", ylim=(0, 1), **kwargs)


def hopPlot(G, **kwargs):
	""" Prints the hop-plot"""
	#hop-plot
	if G.isDirected():
		cc = components.StronglyConnectedComponents(G)
	else:
		cc = components.ConnectedComponents(G)
	cc.run()
	if cc.numberOfComponents() == 1:
		hopPlot = distance.EffectiveDiameter.hopPlot(G, maxDistance=0, k=64, r=7)
	else:
		hopPlot = {}
	plt.xlabel('distance')
	plt.ylabel('fraction of connected nodes')
	plt.ylim([0,1.02])
	plt.plot(list(hopPlot.keys()), list(hopPlot.values()), **kwargs)
	#plt.show()
