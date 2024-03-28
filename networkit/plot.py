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

def nodeAttributes(G, attribute=None):
	""" 
	nodeAttributes(G, sorted=True)
	
	General plotting function for a node attributes using matplotlib.
	
	Parameters
	----------
	G : `*kargs`
		Input Graph of which node attributes are being plotted. 
	attribute : nk.graph.NodeAttribute or tuple(nk.graph.NodeAttribute)
		(tuple of) nk.graph.NodeAttribute attached to the graph.
	"""
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	if attribute is None:
		raise "Error, call with no attribute given." 

	if type(attribute) == graph.NodeAttribute:
		# single attribute, plots number of nodes for each attribute value
		attrFreq= {} 
		for i in G.iterNodes():
			try:
				if attribute[i] not in attrFreq:
					attrFreq[attribute[i]] = 1
				else:
					attrFreq[attribute[i]] += 1
			except ValueError:
				continue
		x,y = zip(*sorted(attrFreq.items()))
		fig, ax = plt.subplots()	
		ax.bar(x, y)
		ax.title.set_text("1-Dim Node Attribute Distribution")
		ax.set_xticks(x)
		ax.set_yticks(y)
		title = str(attribute.getName().decode())
		ax.set_xlabel(title)
		ax.set_ylabel("Number of Nodes")
		plt.show()

	elif type(attribute) == tuple:
		# two attributes, plots 2-Dimensional distribution of nodes
		attributeList = list(attribute)
		plotList = []
		for attr in attributeList:
			data=[0]*G.numberOfNodes()
			for i in G.iterNodes():
				try:
					data[i]=attr[i]
				except ValueError:
					continue
			plotList.append(data)	
		plt.plot(plotList[0], plotList[1], "ro")
		plt.title("2-Dim Node Attribute Distribution")
		titleX = str(attribute[0].getName().decode())
		titleY = str(attribute[1].getName().decode())
		plt.xlabel(titleX)
		plt.ylabel(titleY)
		plt.show()		
	else:
		raise "Error, attribute has wrong type, call with nk.graph.NodeAttribute or tuple(nk.graph.NodeAttribute)." 

def degreeDistribution(G, *args, **kwargs):
	"""
	degreeDistribution(G, *args, **kwargs)

	Plots the degree distribution of the given network using matplotlib.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`*args` : list()
		Additional *args parameter passed to matplotlib.pyplot.bar.
	`**kwargs` : dict()
		Additional **kwargs parameter passed to matplotlib.pyplot.bar
	"""
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	dd = [0] * (graphtools.maxDegree(G) + 1)
	nodes = []
	for i in range(G.numberOfNodes()):
		dd[G.degree(i)] += 1
	for i in range(len(dd)):
		nodes.append(i)

	plt.bar(nodes, dd, *args, **kwargs)
	plt.title("Degree Distribution")
	plt.xlabel("Degree")
	plt.ylabel("Number of Nodes")
	plt.show()

def connectedComponentsSizes(G, relativeSizes=True):
	""" 
	connectedComponentsSizes(G, relativeSizes=True)
	
	Plot the size distribution of connected components as a pie chart using matplotlib. 
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	relativeSizes : bool, optional
		If relativeSizes is set to True, the component sizes in the pie chart will
		correlate with their number of nodes. Default: True
	"""
	if not have_seaborn:
		raise MissingDependencyError("seaborn")
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	csizes = components.ConnectedComponents(G).run().getComponentSizes()
	data = sorted(list(csizes.values()), reverse=True)
	colors = seaborn.color_palette("Set2", 10)
	total = sum(data)
	# explode the largest component pie piece
	maxi = data.index(max(data))
	explode = [0 for i in range(len(data))]
	explode[maxi] = 0.1
	# plot
	plt.figure(figsize=(5,5))
	plt.title("Size of Connected Components")

	def filter_autopct(pct):
		return ('%1.f%%' % pct) if pct > 5 else ''

	if relativeSizes: 
		plt.pie(data, colors=colors, autopct=filter_autopct, explode=explode)
	else: 
		plt.pie(data, colors=colors, autopct=lambda p: '{:.0f}'.format(p * total / 100), explode=explode)

def coreDecompositionSequence(G, *args, **kwargs):
	""" 
	coreDecompositionSequence(G, *args, **kwargs)
	
	Plots the core decomposition sequence of G, i.e. the size of the k-shell for the core number k using matplotlib.
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`*args` : list()
		Additional *args parameter passed to matplotlib.pyplot.bar.
	`**kwargs` : dict()
		Additional **kwargs parameter passed to matplotlib.pyplot.bar
	"""
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	shells = centrality.CoreDecomposition(G).run().getPartition().subsetSizes()
	k=[]
	for i in range(len(shells)):
		k.append(i+1)

	plt.bar(k, shells, *args, **kwargs)
	plt.title("Size of Core Decomposition K-Shells")
	plt.xticks(k)
	plt.yticks(shells)
	plt.xlabel("K-core decomposition(k)")
	plt.ylabel("Size of k-shell")
	plt.show()

def clusteringPerDegree(G):
	""" 
	clusteringPerDegree(G, **kwargs)
	
	Plots the local clustering coefficient for all degrees that exist in the given graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	if not have_seaborn:
		raise MissingDependencyError("seaborn")
	if not have_pandas:
		raise MissingDependencyError("pandas")
	degs = centrality.DegreeCentrality(G).run().scores()
	cc = centrality.LocalClusteringCoefficient(G).run().scores()
	data = pandas.DataFrame({"deg": degs, "cc" : cc})
	data = data.groupby("deg", as_index=False).mean()
	jointplot = seaborn.jointplot(data, x="deg", y="cc",kind="reg", ylim=(0, 1))

def hopPlot(G, *args, **kwargs):
	""" 
	hopPlot(G, **kwargs)
	
	Prints the hop-plot using matplotlib.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`*args` : list()
		Additional *args parameter passed to matplotlib.pyplot.bar.
	`**kwargs` : dict()
		Additional **kwargs parameter passed to matplotlib.pyplot.bar
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
		hopPlot = distance.HopPlotApproximation(G).run().getHopPlot()
	else:
		hopPlot = {}
	distances = list(hopPlot.values())	
	for i in range(len(distances)):
		distances[i]*=100 #turn fractions into percentages
	plt.title("Hop Plot Approximation")
	plt.xlabel('Distance (d)')
	plt.ylabel('Percentage of connected nodes (g(d))')
	plt.ylim([0,102])
	plt.plot(list(hopPlot.keys()), distances, *args, **kwargs)
	plt.show()
