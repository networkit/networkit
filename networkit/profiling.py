#
# file: profiling.py
# author: Mark Erb
#

from networkit import *

import multiprocessing
from IPython.core.display import *
from urllib.parse import quote
import io
import sys

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
			
			
try:
	__IPYTHON__
except:
	raise ImportError("module has been loaded outside of \"IPython\"")

	
def readfile(postfix):
	with open(__file__[:__file__.rfind(".py")] + "." + postfix, "r") as file:
		return " ".join(file.read().split())
			

def __initHeader(tag, type, data):
	result = """
		{
			var element = document.getElementById('NetworKit_""" + tag + """');
			if (element) {
				element.parentNode.removeChild(element);
			}
			element = document.createElement('""" + tag + """');
			element.type = 'text/""" + type + """';
			element.innerHTML = '""" + data + """';
			element.setAttribute('id', 'NetworKit_""" + tag + """');
			document.head.appendChild(element);
		}
	"""
	return result
		
		
def __initOverlay(name, data):
	result = """
		{
			var element = document.getElementById('NetworKit_""" + name + """');
			if (element) {
				element.parentNode.removeChild(element);
			}
			element = document.createElement('div');
			element.innerHTML = '<div id="NetworKit_""" + name + """_Toolbar_Top"><div class="button icon-close" id="NetworKit_""" + name + """_Close" /></div>""" + data + """';
			element.setAttribute('id', 'NetworKit_""" + name + """');
			document.body.appendChild(element);
			document.getElementById('NetworKit_""" + name + """_Close').onclick = function (e) {
				document.getElementById('NetworKit_""" + name + """').style.display = 'none';
			}
		}
	"""	
	return result


display_html(
	HTML("""
		<script type="text/javascript">
		<!--
			""" + __initHeader("script", "javascript", readfile("js"))  + """
			""" + __initHeader("style",  "css",        readfile("css")) + """
			""" + __initOverlay("Overlay", readfile("overlay.html")) + """
		-->
		</script>
	""")
)


class Worker(multiprocessing.Process):
	def __init__(self, tasks, results):
		multiprocessing.Process.__init__(self)
		self.__tasks = tasks
		self.__results = results
		
	def run(self):
		while True:
			task = self.__tasks.get()
			if task is None:
				self.__tasks.task_done()
				break
			result = task.run()
			self.__tasks.task_done()
			self.__results.put(result)
		
	
class Plot_Task(object):
	def __init__(self, name, index, data):
		self.__name = name
		self.__index = index
		self.__data = data
	
	
	def run(self):
		plt.ioff()	
		xLabel = "x-Axis"
		yLabel = "y-Axis"
		fig = plt.figure()
		n, bins, patches = plt.hist(self.__data, 50, normed=1, facecolor='green', alpha=0.75)
		plt.xlabel(xLabel)
		plt.ylabel(yLabel)
		plt.grid(True)
		fig.tight_layout()
		imgdata = io.StringIO()
		fig.savefig(imgdata, format='svg')
		plt.close(fig)
		plaintext = imgdata.getvalue()
		plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
		encoded = quote(plaintext, safe='');
		return (self.__name, self.__index, encoded)
		
	
class Profiling:
	__TOKEN = object();		
	__pageCount = 0
	__verbose = False
	__parallel = multiprocessing.cpu_count() * 2
		
		
	def __init__(self, G, token):
		if token is not self.__TOKEN:
			raise ValueError("call create(G) to create an instance")
		self.__G = G
		self.__measures = {}
		
	
	@classmethod
	def create(cls, G):
		result = cls(G, cls.__TOKEN)
		
		for parameter in [ 
			(centrality.DegreeCentrality, 			(G, )),
			(centrality.CoreDecomposition, 			(G, )),
			(centrality.LocalClusteringCoefficient, (G, )),
			(centrality.PageRank, 					(G, )),
			(centrality.KPathCentrality,			(G, )),
			(centrality.KatzCentrality,				(G, )),
			(centrality.ApproxBetweenness2,			(G, max(42, G.numberOfNodes() / 1000), False))
		]: result.__addMeasure(parameter)
		
		result.__loadMeasures()
		return result;
	

	@classmethod
	def setVerbose(cls, verbose):
		cls.__verbose = verbose
	
	
	@classmethod
	def getVerbose(cls):
		return cls.__verbose
	
	
	@classmethod
	def setParallel(cls, parallel):
		if (parallel < 1):
			raise ValueError("parallel < 1");
		cls.__parallel = parallel
	

	@classmethod
	def getParallel(cls):
		return cls.__parallel


	def show(self):
		pageIndex = self.__pageCount
		tasks = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		
		if self.__verbose:
			timerAll = stopwatch.Timer()
		
		workers = [Worker(tasks, results) for i in range(self.__parallel)]
		for w in workers:
			w.deamon = True
			w.start()
			
		numberOfTasks = 0
		for key in self.__measures:
			tasks.put(Plot_Task(key, 0, self.__measures[key]["data"]))
			numberOfTasks += 1
		for i in range(self.__parallel):
			tasks.put(None)
		tasks.join()
		
		while(numberOfTasks):
			(name, index, image) = results.get()
			self.__measures[name]["image"] = image
			numberOfTasks -= 1
		
		centralities = ""
		for key in self.__measures:
			image = self.__measures[key]["image"]
			centralities = centralities + "<div class=\"Plot\" title=\"" + key + "\" data-image=\"data:image/svg+xml;utf8," + image + "\" />"
			
		result = readfile("profile.html")
		result = result.format(**locals());
		display_html(HTML(result))
		
		if self.__verbose:
			print("\ntotal time: {:.2F} s".format(timerAll.elapsed))
	
		self.__pageCount = self.__pageCount + 1
	
	
	def __addMeasure(self, args):
		(measureClass, parameters) = args
		measureName = measureClass.__name__
		measure = {}
		measure["index"] = len(self.__measures)
		measure["class"] = measureClass
		measure["parameters"] = parameters
		self.__measures[measureName] = measure
	
	
	def __loadMeasures(self):
		if self.__verbose:
			timerAll = stopwatch.Timer()
		
		for name in self.__measures:
			measure = self.__measures[name]
			instance = measure["class"](*measure["parameters"])
			if self.__verbose:
				print(name + ": ", end="", flush=True)
			timerInstance = stopwatch.Timer()
			instance.run()
			elapsed = timerInstance.elapsed
			if self.__verbose:
				print("{:.2F} s".format(elapsed))
			measure["data"] = instance.scores()
			measure["time"] = elapsed
		if self.__verbose:
			print("\ntotal time: {:.2F} s".format(timerAll.elapsed))	
			
	
# class Plot:
	# __metaclass__ = ABCMeta
			
	# def init(self, xScale, yScale):
		# plt.clf();
		# fig, ax = plt.subplots()
		# if xScale:
			# ax.set_xscale('log')
		# if yScale:
			# ax.set_yscale('log')


	# def createImageURI(self):
		# fig = plt.gcf()
		# fig.tight_layout()
		# imgdata = io.StringIO()
		# fig.savefig(imgdata, format='svg')
		# plaintext = imgdata.getvalue()
		# plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
		# encoded = quote(plaintext, safe='');
		# encoded = "data:image/svg+xml;utf8," + encoded; 
		# return encoded;
		
			
	# def plotSet(self, data, xLabel, yLabel):
		# result = ""
		# for xScale in range(0, 2):
			# for yScale in range(0, 2):
				# result += self.plot(data, xLabel, yLabel, xScale, yScale)
				# if (not(xScale and yScale)):
					# result += "|"
		# return result;
			
		
	# @abstractmethod 
	# def plot(self, data, xLabel, yLabel, xScale, yScale): pass
			
			
# class Histogram(Plot):
	# def plot(self, data, xLabel, yLabel, xScale, yScale):
		# self.init(xScale, yScale)
		# n, bins, patches = plt.hist(data, 50, normed=1, facecolor='green', alpha=0.75)
		# plt.xlabel(xLabel)
		# plt.ylabel(yLabel)
		# plt.grid(True)
		# return self.createImageURI()
				
				
# class HistogramSeaborn(Plot):
	# def plot(self, data, xLabel, yLabel, xScale, yScale):
		# self.init(xScale, yScale)
		# sns.distplot(data);
		# plt.xlabel(xLabel)
		# plt.ylabel(yLabel)
		# plt.grid(True)
		# return self.createImageURI()
	
# def asImage(plotFunction, plotArgs=[], plotKwargs={}, size=(8,6)):
	# """
	# Call any plot function with the given argument and return the image in an HTML <img> tag.
	# """
	# plt.figure(figsize=size)
	# plotFunction(*plotArgs, **plotKwargs)
	# # Get a handle for the plot that was just generated
	# fig = Gcf.get_all_fig_managers()[-1].canvas.figure
	# # Generate a data URL for the image
	# imageData = "data:image/png;base64,{0}".format(b64encode(print_figure(fig)).decode("utf-8"))
	# # Remove the plot from the list of plots for the current cell
	# Gcf.destroy_fig(fig)
	# # generate img tag
	# image = "<img src='{0}'>".format(imageData)
	# return image

# def computeNetworkProperties(G):
	# """
	# """
	# networkProperties = [
			# ["nodes, edges", "{0}, {1}".format(G.numberOfNodes(), G.numberOfEdges())],
			# ["directed?", "{0}".format(G.isDirected())],
			# ["weighted?", "{0}".format(G.isWeighted())],
			# #["density", "{0:.6f}".format(properties.density(G))],
			# ["diameter range", "{0}".format(properties.Diameter.estimatedDiameterRange(G, error=0.1))]
		# ]
	# return networkProperties


# def computeNodePartitions(G):
	# components = properties.components(G)
	# communities = community.detectCommunities(G)

	
# def computeNodeProperties(G):
	# # degree
	# degree = properties.degreeSequence(G)
	# # coreness
	# core = centrality.CoreDecomposition(G).run().scores()
	# # local clustering coefficient
	# clustering = centrality.LocalClusteringCoefficient(G).run().scores()
	# # betweenness
	# nSamples = max(42, G.numberOfNodes() / 1000)
	# betweenness = centrality.ApproxBetweenness2(G, nSamples, normalized=True).run().scores()
	# # pagerank
	# pagerank = centrality.PageRank(G).run().scores()
	# # k-Path centrality
	# kpath = centrality.KPathCentrality(G).run().scores()
	# # Katz centrality
	# katz = centrality.KatzCentrality(G).run().scores()
	# # package node properties in DataFrame
	# nodeProperties = pandas.DataFrame({"degree": degree,
	 									# "core": core,
										# "clustering": clustering,
										# "betweenness": betweenness,
										# "pagerank": pagerank,
										# "kpath": kpath,
										# "katz": katz})
	# return nodeProperties


# def computeNodePropertyCorrelations(nodeProperties, method="spearman"):
	# return nodeProperties.corr(method=method)


# def plotNodePropertyCorrelations(nodeProperties, figsize=(8,8), method="spearman"):
    # cmap = seaborn.diverging_palette(220, 20, as_cmap=True)
    # f, ax = plt.subplots(figsize=figsize)
    # print("correlating"); sys.stdout.flush()
    # seaborn.corrplot(nodeProperties, cmap=cmap, method=method)
    # f.tight_layout()



# def profile(G):
	# """
	# Output profile page of network as HTML
	# """
	# global loaded
	# # settings
	# defaultSize = (5,2)
	# histArgs = {"bins" : 100, "figsize" : (12,8)}

	# # compute global network attributes
	# networkProperties = computeNetworkProperties(G)
	# networkPropertiesTable = tabulate.tabulate(networkProperties, tablefmt="html")

	# hopPlot = asImage(plot.hopPlot, plotArgs=[G], size=defaultSize)

	# # compute node properties
	# nodeProperties = computeNodeProperties(G)


	# # compute figures
	# (plaw, _, gamma) = properties.degreePowerLaw(G)

	# # compute images
	# ddPlot = asImage(plot.degreeDistribution, plotArgs=[G], size=defaultSize)
	# ddHist = asImage(nodeProperties["degree"].hist, plotKwargs=histArgs, size=defaultSize)
	# dd = asSlideShow(ddPlot, ddHist)
				  
	# ccPlot = asImage(plot.nodeProperty, plotKwargs={"data" : nodeProperties["clustering"], "label": "local clustering coefficient"}, size=defaultSize)
	# ccHist = asImage(nodeProperties["clustering"].hist, plotKwargs=histArgs, size=defaultSize)
	# compPlot = asImage(plot.connectedComponentsSizes, [G], size=(1,1))

	# html = readfile("html")

	# result = html.format(**locals())
	# return HTML(result)