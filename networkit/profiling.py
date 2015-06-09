from networkit import *
import tabulate
import pandas
import seaborn
import matplotlib.pyplot as plt
from matplotlib._pylab_helpers import Gcf
from IPython.core.pylabtools import print_figure
from base64 import b64encode
from IPython.core.display import HTML


profileTemplate = """
	<style media="screen" type="text/css">
	#wrapper {
	    width: 500px;
	    border: 1px solid black;
	}
	#first {
	    width: 300px;
	    border: 1px solid red;
	}
	#second {
	    border: 1px solid green;
	}
	</style>

	<div id="page">
	<h1>Network Profile</h1>

	<h2>Network Properties</h2>
		<div id="wrapper">
			<div id="first">
				{networkPropertiesTable}
			</div>
			<div id="second">
				{hopPlot}
			</div>
		</div>

	<h2>Network Partitions</h2>

	<h2>Node Centrality Measures</h2>

	<h3>Degree</h3>
	{ddPlot}
	{ddHist}

	<h3>Local Clustering Coefficient</h3>
	{ccPlot}
	{ccHist}




	power law distribution: {plaw} {gamma}

	{compPlot}
	</div>
"""


def asImage(plotFunction, plotArgs=[], plotKwargs={}, size=(8,6)):
	"""
	Call any plot function with the given argument and return the image in an HTML <img> tag.
	"""
	plt.figure(figsize=size)
	plotFunction(*plotArgs, **plotKwargs)
	# Get a handle for the plot that was just generated
	fig = Gcf.get_all_fig_managers()[-1].canvas.figure
	# Generate a data URL for the image
	imageData = "data:image/png;base64,{0}".format(b64encode(print_figure(fig)).decode("utf-8"))
	# Remove the plot from the list of plots for the current cell
	Gcf.destroy_fig(fig)
	# generate img tag
	image = "<img src='{0}'\>".format(imageData)
	return image

def computeNetworkProperties(G):
	"""
	"""
	networkProperties = [
			["nodes, edges", "{0}, {1}".format(G.numberOfNodes(), G.numberOfEdges())],
			["directed?", "{0}".format(G.isDirected())],
			["weighted?", "{0}".format(G.isWeighted())],
			["density", "{0:.6f}".format(properties.density(G))],
			["diameter range", "{0}".format(properties.Diameter.estimatedDiameterRange(G, error=0.1))]
		]
	return networkProperties


def computeNodePartitions(G):
	components = properties.components(G)
	communities = community.detectCommunities(G)

def computeNodeProperties(G):
	# degree
	degree = properties.degreeSequence(G)
	# coreness
	core = centrality.CoreDecomposition(G).run().scores()
	# local clustering coefficient
	clustering = centrality.LocalClusteringCoefficient(G).run().scores()
	# betweenness
	nSamples = max(42, G.numberOfNodes() / 1000)
	betweenness = centrality.ApproxBetweenness2(G, nSamples, normalized=True).run().scores()
	# pagerank
	pagerank = centrality.PageRank(G).run().scores()
	# k-Path centrality
	kpath = centrality.KPathCentrality(G).run().scores()
	# Katz centrality
	katz = centrality.KatzCentrality(G).run().scores()
	# package node properties in DataFrame
	nodeProperties = pandas.DataFrame({"degree": degree,
	 									"core": core,
										"clustering": clustering,
										"betweenness": betweenness,
										"pagerank": pagerank,
										"kpath": kpath,
										"katz": katz})
	return nodeProperties


def computeNodePropertyCorrelations(nodeProperties, method="spearman"):
	return nodeProperties.corr(method=method)


def plotNodePropertyCorrelations(nodeProperties, figsize=(8,8), method="spearman"):
    cmap = seaborn.diverging_palette(220, 20, as_cmap=True)
    f, ax = plt.subplots(figsize=figsize)
    print("correlating"); sys.stdout.flush()
    seaborn.corrplot(nodeProperties, cmap=cmap, method=method)
    f.tight_layout()



def profile(G):
	"""
	Output profile page of network as HTML
	"""
	# settings
	defaultSize = (5,2)
	histArgs = {"bins" : 100, "figsize" : (12,8)}

	# compute global network attributes
	networkProperties = computeNetworkProperties(G)
	networkPropertiesTable = tabulate.tabulate(networkProperties, tablefmt="html")

	hopPlot = asImage(plot.hopPlot, plotArgs=[G], size=defaultSize)

	# compute node properties
	nodeProperties = computeNodeProperties(G)


	# compute figures
	(plaw, _, gamma) = properties.degreePowerLaw(G)

	# compute images
	ddPlot  = asImage(plot.degreeDistribution, plotArgs=[G], size=defaultSize)
	ddHist = asImage(nodeProperties["degree"].hist, plotKwargs=histArgs, size=defaultSize)
	ccPlot = asImage(plot.nodeProperty, plotKwargs={"data" : nodeProperties["clustering"], "label": "local clustering coefficient"}, size=defaultSize)
	ccHist = asImage(nodeProperties["clustering"].hist, plotKwargs=histArgs, size=defaultSize)
	compPlot = asImage(plot.connectedComponentsSizes, [G], size=(1,1))

	page = HTML(profileTemplate.format(**locals()))
	return page
