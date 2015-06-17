from networkit import *
import tabulate
import pandas
import seaborn
import powerlaw
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
	partitions = {	"components":	properties.components(G),
					"communities":	community.detectCommunities(G)}
	# TODO: shells
	for (partitionName, partition) in partitions.items():
		partition.subsetSizes()

	raise NotImplementedError("TODO")

def computeNodeCentralities(G):
	(n, m) = G.size()

	# TODO: normalization?

	nodeCentralityAlgos = {
							"degree":		(centrality.DegreeCentrality, 			(G, )),
							"coreness":		(centrality.CoreDecomposition, 			(G, )),
							"clustering":	(centrality.LocalClusteringCoefficient, (G, )),
							"pagerank":		(centrality.PageRank, 					(G, )),
							"kpath":		(centrality.KPathCentrality,			(G, )),
							"katz":			(centrality.KatzCentrality,				(G, )),
							"betweenness":	(centrality.ApproxBetweenness2,			(G, max(42, n / 1000))),
							"closeness":	(centrality.Closeness,					(G, max(42, n / 1000)))
							}

	centralityScores = {}
	for (algoName, (algoClass, params)) in nodeCentralityAlgos.items():
		algo = algoClass(*params)
		t = stopwatch.Timer()
		algo.run()
		print(algoName, ": ", "{:.2E}".format(t.elapsed))
		centralityScores[algoName] = algo.scores()
	nodeCentralities = pandas.DataFrame(centralityScores)
	return nodeCentralities

def powerLawStats(centralities):
	powerLawStats = {}
	for (centralityName, centralityScores) in centralities.items():
		fit = powerlaw.Fit(centralityScores)
		R, p = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)
		gamma = fit.alpha
		powerLawStats[centralityName] = ((R > 0), R, gamma)
	return powerLawStats

def computeRankCorrelations(centralities : pandas.DataFrame, method="spearman"):
	return centralities.corr(method=method)


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
	nodeProperties = perties(G)


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


class Profile:
	""" This class computes and presents a structural profile of a networks"""

	def __init__(self, G, settings={}):
		self.G = G
		self.settings = settings

	def compute(self):
		self.nodeCentralities = computeNodeCentralities(self.G)
		raise NotImplementedError("TODO")

	def getPage(self):
		raise NotImplementedError("TODO")

	def getAttributeVector(self):
		raise NotImplementedError("TODO:")
