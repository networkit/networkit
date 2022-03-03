from networkit import *
import numpy as np
import operator
import warnings
from .support import MissingDependencyError
try:
	import matplotlib.pyplot as plt
except ImportError:
	have_plt = False
else:
	have_plt = True
try:
	import pandas
except ImportError:
	have_pandas = False
else:
	have_pandas = True
try:
	import seaborn
	seaborn.set_style("whitegrid")
except ImportError:
	have_seaborn = False
else:
	have_seaborn = True

def nodeProperty(data, label, sorted=True, yscale="linear", xscale="linear"):
	""" 
	nodeProperty(data, label, sorted=True, yscale="linear", xscale="linear")
	
	General plotting function for a node property using matplotlib.
	
	Parameters
	----------
	data : `*kargs`
		Input data for matplotlib.plot function. Refer to offical documentation for details: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
	label : list(str)
		Label for data points.
	sorted : bool, optional
		Indicates whether the plotted data points should be sorted. Default: True
	yscale : str, optional
		Indicates the scaling of y-axis. Default: True
	xscale : str, optional
		Indicates the scaling of x-axis. Default: True
	"""
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	plt.yscale(yscale)
	plt.xscale(xscale)
	plt.xlabel("nodes")
	plt.ylabel(label)
	if sorted:
		data.sort(reverse=True)
	plt.plot(data)


def degreeDistribution(G, **kwargs):
	"""
	degreeDistribution(G, **kwargs)

	Plots the degree distribution of the given network using matplotlib.
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`**kwargs` : `**kwargs`
		Input data currently not used.
	"""
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	dd = properties.degreeDistribution(G)
	plt.yscale("symlog")
	plt.xscale("log")
	plt.xlabel("nodes")
	plt.ylabel("degree")
	plt.plot(dd)

def connectedComponentsSizes(G, **kwargs):
	""" 
	connectedComponentsSizes(G, **kwargs)
	
	Plot the size distribution of connected components as a pie chart using matplotlib. 
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`**kwargs` : `**kwargs`
		Input parameter currently not used.
	"""
	if not have_seaborn:
		raise MissingDependencyError("seaborn")
	if not have_plt:
		raise MissingDependencyError("matplotlib")
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
	""" 
	coreDecompositionSequence(G, **kwargs)
	
	Plots the core decomposition sequence of G, i.e. the size of the k-shell for the core number k using matplotlib.
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`**kwargs` : `**kwargs`
		Input parameter currently not used.
	"""
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	if not have_pandas:
		raise MissingDependencyError("pandas")
	shells = centrality.CoreDecomposition(G).run().shells()
	data = pandas.DataFrame({"k": range(len(shells)), "n_k": [len(shell) for shell in shells]})
	plt.xlabel("core number")
	plt.ylabel("number of nodes")
	plt.plot(data["k"], data["n_k"], **kwargs)


def clusteringPerDegree(G, **kwargs):
	""" 
	clusteringPerDegree(G, **kwargs)
	
	Plots the local clustering coefficient for nodes with specific degree.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`**kwargs` : `**kwargs`
		Input parameter currently not used.
	"""
	if not have_seaborn:
		raise MissingDependencyError("seaborn")
	if not have_pandas:
		raise MissingDependencyError("pandas")
	degs = centrality.DegreeCentrality(G).run().scores()
	cc = centrality.LocalClusteringCoefficient(G).run().scores()
	data = pandas.DataFrame({"deg": degs, "cc" : cc})
	data = data.groupby("deg", as_index=False).mean()
	jointplot = seaborn.jointplot("deg", "cc", data, kind="reg", ylim=(0, 1), **kwargs)


def hopPlot(G, **kwargs):
	""" 
	hopPlot(G, **kwargs)
	
	Prints the hop-plot using matplotlib.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`**kwargs` : `**kwargs`
		Input parameter currently not used.
	"""
	#hop-plot
	if not have_plt:
		raise MissingDependencyError("matplotlib")
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
