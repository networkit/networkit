#
# file: profiling.py
# author: Mark Erb
#

from networkit import *

import multiprocessing
from IPython.core.display import *
from urllib.parse import quote
import io
import math

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



def hoelderMean(sample, n, p):
	result = 0
	for i in range(n):
		result += sample[i] ** p
	result /= n
	result **= 1 / p
	return result
	
	
def momentum(sample, n, arithmeticMean, s_n, p):
	result = 0
	for i in range(n):
		result += ((sample[i] - arithmeticMean) / s_n) ** p
	result /= n
	return result


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
			result = (task.getType(), task.getName(), task.run())
			self.__tasks.task_done()
			self.__results.put(result)

			
class Stat_Task(object):
	def __init__(self, name, params):
		self.__name = name
		self.__params = params
	
	def getName(self):
		return self.__name
		
	def getType(self):
		return "Stat"
	
	def run(self):
		sample = self.__params
		n = len(sample)
		sampleSorted = sorted(sample)
	
		results = {}
		results["Location"] =  {}
		results["Dispersion"] =  {}
		results["Shape"] =  {}
		results["Binning"] = {}
		results["Size"] = n
		
		results["Location"]["Arithmetic Mean"] = arithmeticMean = hoelderMean(sample, n, 1)
		results["Location"]["Quadratic Mean"] = quadraticMean = hoelderMean(sample, n, 2)
		
		def funcVariance():
			result = 0
			for i in range(n):
				result += (sample[i] - arithmeticMean) ** 2
			result /= n - 1
			return result
		results["Dispersion"]["Variance"] = variance = funcVariance()
		
		def funcStandardDeviation():
			result = math.sqrt(variance)
			return result
		results["Dispersion"]["Standard Deviation"] = s_n = funcStandardDeviation()
		
		def funcCoefficientOfVariation():
			result = s_n / arithmeticMean
			return result
		results["Dispersion"]["Coefficient Of Variation"] = c_v = funcCoefficientOfVariation()
		
		def funcMin():
			result = sampleSorted[0]
			return result
		results["Location"]["Min"] = min = funcMin()
			
		def funcMax():
			result = sampleSorted[n-1]
			return result	
		results["Location"]["Max"] = max = funcMax()
		
		def funcAlphaQuartile(alpha):
			k_real = (alpha * n)
			k = math.floor(k_real)
			if (k != k_real) or (k < 1):
				result = sampleSorted[(k-1)+1]
			else:
				result = 0.5 * (sampleSorted[(k-1)] + sampleSorted[(k-1)+1])
			return result
		results["Location"]["1st Quartile"] = Q1 = funcAlphaQuartile(0.25)
		results["Location"]["Median"] = median = funcAlphaQuartile(0.5)
		results["Location"]["3rd Quartile"] = Q3 = funcAlphaQuartile(0.75)
		
		def funcAlphaTrimmedMean(alpha):
			k = math.floor(alpha * n)
			i = k+1
			result = 0
			while (i < n-k+1):
				result += sampleSorted[(i-1)]
				i += 1
			result /= n	- 2*(k)
			return result
		results["Location"]["Interquartile Mean"] = IQM = funcAlphaTrimmedMean(0.25)
		
		def funcIQR():
			result = Q3 - Q1
			return result
		results["Dispersion"]["Interquartile Range"] = IQR = funcIQR()
		
		def funcSampleRange():
			result = max - min
			return result
		results["Dispersion"]["Sample Range"] = sampleRange = funcSampleRange()
		
		def funcSkewnessYP():
			result = 3 * (arithmeticMean - median) / s_n
			return result
		results["Shape"]["Skewness YP"] = skewness_yp = funcSkewnessYP()
		
		def funcSkewnessM():
			result = momentum(sample, n, arithmeticMean, s_n, 3)
			return result
		results["Shape"]["Skewness M"] = skewnewss_m = funcSkewnessM()
		
		def funcKurtosis():
			result = momentum(sample, n, arithmeticMean, s_n, 4) - 3
			return result
		results["Shape"]["Kurtosis"] = kurtosis = funcKurtosis()
		
		def funcNumberOfBins():
			result = math.sqrt(n)
			if (result < 5):
				result = 5
			elif (result > 20):
				result = 20
			return result
		results["Binning"]["Number"] = k_Bins = funcNumberOfBins()
		
		def funcIntervals():
			result = []
			w = sampleRange / k_Bins
			result.append(min)
			for i in range(1, k_Bins):
				result.append(min + w * i)
			result.append(max)
			return result
		results["Binning"]["Intervals"] = intervals  = funcIntervals()
		
		def funcBinAbsoluteFrequencies():
			result = [0]
			index = 0
			for i in range(n):
				value = sampleSorted[i]
				if intervals[index + 1] < value:
					result.append(0)
					index += 1
				result[index] += 1
			return result
		results["Binning"]["Absolute Frequencies"] = absoluteFrequencies = funcBinAbsoluteFrequencies()
		
		
		def funcBinRelativeFrequencies():
			result = []
			for H in absoluteFrequencies:
				result.append(H / n)
			return result
		results["Binning"]["Relative Frequencies"] = relativeFrequencies = funcBinRelativeFrequencies()
		
		def funcMode():
			index = 0
			max = 0
			for i in range(len(absoluteFrequencies)):
				if absoluteFrequencies[i] > max:
					max = absoluteFrequencies[i]
					index = i
			result = (intervals[index]+intervals[index+1]) / 2
			return result
		results["Binning"]["Mode"] = mode = funcMode()
		
		def funcSkewnessKPM():
			result = (arithmeticMean - mode) / s_n
			return result
		results["Shape"]["Skewness KPM"] = skewnewss_kpm = funcSkewnessKPM()
		
		return results

		
class Corr_Task(object):
	def __init__(self, name, params):
		self.__name = name
		self.__params = params
	
	def getName(self):
		return self.__name
		
	def getType(self):
		return "Stat"
	
	def run(self):
		(sample_1, sample_2) = self.__params
		n = len(sample_1)
		assert (n == len(sample_2)), "sample size are not equal"
		
		results = {}
		
		def funcRanged():
			result = []
			for i in range(n):
				result.append([sample_1[i], -1, sample_2[i], -1])
			for i in [0, 2]:
				result.sort(key=lambda x: x[i])
				value = result[0][i]
				sum = 0
				length = 0
				for j in range(n):
					if value == result[j][i]:
						sum += (j+1)
						length += 1
					else:
						sum /= length
						for k in range(length):
							result[j-k-1][i+1] = sum
						value = result[j][i]
						sum = (j+1)
						length = 1
				sum /= length
				for j in range(length):
					result[n-j-1][i+1] = sum
			return result
		ranged = funcRanged()
		
		def funcRangedArithmeticMean(p):
			result = 0
			for i in range(n):
				result += ranged[i][p]
			result /= n
			return result
		rangedArithmeticMean_1 = funcRangedArithmeticMean(1)
		rangedArithmeticMean_2 = funcRangedArithmeticMean(3)
		
		def funcRangedCovariance():
			result = 0
			for i in range(n):
				result += (ranged[i][1] - rangedArithmeticMean_1) * (ranged[i][3] - rangedArithmeticMean_2)
			result /= n
			return result
		rangedCovariance = funcRangedCovariance()
		
		def funcRangedVariance(p, rangedArithmeticMean):
			result = 0
			for i in range(n):
				result += (ranged[i][p] - rangedArithmeticMean) ** 2
			result /= n
			return result
		rangedVariance_1 = funcRangedVariance(1, rangedArithmeticMean_1)
		rangedVariance_2 = funcRangedVariance(3, rangedArithmeticMean_2)
		
		def funcRangedStandardDeviation(rangedVariance):
			result = math.sqrt(rangedVariance)
			return result
		rangedStandardDeviation_1 = funcRangedStandardDeviation(rangedVariance_1)
		rangedStandardDeviation_2 = funcRangedStandardDeviation(rangedVariance_2)
		
		def funcSpearmansCorrelationCoefficient():
			result = rangedCovariance / (rangedStandardDeviation_1 * rangedStandardDeviation_2)
			return result
		results["Spearman's Correlation Coefficient"] = funcSpearmansCorrelationCoefficient()
		
		return results

	
class Plot_Task(object):
	def __init__(self, name, params):
		self.__name = name
		self.__params = params
	
	def getName(self):
		return self.__name
		
	def getType(self):
		return "Plot"
	
	def run(self):
		(index, data) = self.__params
		plt.ioff()	
		xLabel = "x-Axis"
		yLabel = "y-Axis"
		fig = plt.figure()
		n, bins, patches = plt.hist(data, 50, normed=1, facecolor='green', alpha=0.75)
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
		return (index, encoded)
		
	
class Profile:
	__TOKEN = object();		
	__pageCount = 0
	__verbose = False
	__parallel = multiprocessing.cpu_count() * 2 - 1
		
		
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
		
		if self.__verbose:
			timerAll = stopwatch.Timer()
		
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
		numberOfTasks = 0
		tasks = multiprocessing.JoinableQueue()
		results = multiprocessing.Queue()
		workers = [Worker(tasks, results) for i in range(self.__parallel)]
		for w in workers:
			w.deamon = True
			w.start()
			
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
			data = instance.scores()
			tasks.put(Plot_Task(name, (0, data)))
			tasks.put(Stat_Task(name, data))
			numberOfTasks += 2
			measure["time"] = elapsed
		
		while(numberOfTasks):
			(type, name, data) = results.get()
			if (type == "Plot"):
				(index, image) = data
				self.__measures[name]["image"] = image
			elif (type == "Stat"):
				self.__measures[name]["stat"] = data
			numberOfTasks -= 1
		
		for i in range(self.__parallel):
			tasks.put(None)
		tasks.join()
		
		if self.__verbose:
			print("\ntotal time (measures + stats + plots): {:.2F} s".format(timerAll.elapsed))	
			
	
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