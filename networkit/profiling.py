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
	
	
def ranged(sample):
	n = len(sample)
	result = []
	for i in range(n):
		result.append([sample[i], i, -1])	
	result.sort(key=lambda x: x[0])
	value = result[0][0]
	sum = 0
	length = 0
	for i in range(n):
		if value == result[i][0]:
			sum += (i+1)
			length += 1
		else:
			sum /= length
			for j in range(length):
				result[i-j-1][2] = sum
			value = result[i][0]
			sum = (i+1)
			length = 1
	sum /= length
	for i in range(length):
		result[n-i-1][2] = sum
	result.sort(key=lambda x: x[1])
	for i in range(n):
		result[i] = result[i][2]
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
		(sample, sampleSorted, sampleRanged) = self.__params
		n = len(sample)
	
		results = {}
		results["Location"] =  {}
		results["Dispersion"] =  {}
		results["Shape"] =  {}
		results["Binning"] = {}
		results["Properties"] = {}
		
		results["Properties"]["Size"] = n
		
		def funcBesselsCorrection():
			result = n / (n-1)
			return result
		results["Properties"]["Bessel's Correction"] = besselsCorrection = funcBesselsCorrection()
		
		def hoelderMean(sample, p):
			result = 0
			for i in range(n):
				result += sample[i] ** p
			result /= n
			result **= 1 / p
			return result			
		results["Location"]["Arithmetic Mean"] = arithmeticMean = hoelderMean(sample, 1)
		results["Location"]["Arithmetic Mean (Rang)"] = arithmeticMean_Rang = hoelderMean(sampleRanged, 1)
		results["Location"]["Quadratic Mean"] = quadraticMean = hoelderMean(sample, 2)
		
		def funcUncorrectedVariance(sample, arithmeticMean):
			result = 0
			for i in range(n):
				result += (sample[i] - arithmeticMean) ** 2
			result /= n
			return result
		results["Dispersion"]["Uncorrected Variance"] = variance_uncorrected = funcUncorrectedVariance(sample, arithmeticMean)
		results["Dispersion"]["Uncorrected Variance (Rang)"] = variance_Rang_uncorrected = funcUncorrectedVariance(sampleRanged, arithmeticMean_Rang)
		
		def funcVariance(variance_uncorrected):
			result = variance_uncorrected * besselsCorrection
			return result
		results["Dispersion"]["Variance"] = variance = funcVariance(variance_uncorrected)
		results["Dispersion"]["Variance (Rang)"] = variance_Rang = funcVariance(variance_Rang_uncorrected)
		
		
		def funcStandardDeviation(variance):
			result = math.sqrt(variance)
			return result
		results["Dispersion"]["Standard Deviation"] = s_n = funcStandardDeviation(variance)
		results["Dispersion"]["Standard Deviation (Rang)"] = s_n_Rang = funcStandardDeviation(variance_Rang)
		results["Dispersion"]["Uncorrected Standard Deviation"] = s_n_uncorrected = funcStandardDeviation(variance_uncorrected)
		results["Dispersion"]["Uncorrected Standard Deviation (Rang)"] = s_n_Rang_uncorrected = funcStandardDeviation(variance_Rang_uncorrected)
		
		
		def funcCoefficientOfVariation(s_n, arithmeticMean):
			result = s_n / arithmeticMean
			return result
		results["Dispersion"]["Coefficient Of Variation"] = c_v = funcCoefficientOfVariation(s_n, arithmeticMean)
		results["Dispersion"]["Coefficient Of Variation (Rang)"] = c_v_Rang = funcCoefficientOfVariation(s_n_Rang, arithmeticMean_Rang)
		results["Dispersion"]["Uncorrected Coefficient Of Variation"] = c_v = funcCoefficientOfVariation(s_n_uncorrected, arithmeticMean)
		results["Dispersion"]["Uncorrected Coefficient Of Variation (Rang)"] = c_v_Rang = funcCoefficientOfVariation(s_n_Rang_uncorrected, arithmeticMean_Rang)
		
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

		def funcMidRange():
			result = sampleRange / 2
			return result
		results["Dispersion"]["Mid-Range"] = midRange = funcMidRange()
		
		def funcSkewnessYP():
			result = 3 * (arithmeticMean - median) / s_n
			return result
		results["Shape"]["Skewness YP"] = skewness_yp = funcSkewnessYP()
		
		def funcMomentum(p):
			result = 0
			for i in range(n):
				result += ((sample[i] - arithmeticMean) / s_n) ** p
			result /= n
			return result
		
		def funcSkewnessM():
			result = funcMomentum(3)
			return result
		results["Shape"]["Skewness M"] = skewnewss_m = funcSkewnessM()
		
		def funcKurtosis():
			result = funcMomentum(4) - 3
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
		(sample_1, sampleRanged_1, stat_1, sample_2, sampleRanged_2, stat_2) = self.__params
		n = len(sample_1)
		assert (n == len(sample_2)), "sample size are not equal"
		
		arithmeticMean_1_Rang = stat_1["Location"]["Arithmetic Mean (Rang)"]
		arithmeticMean_2_Rang = stat_2["Location"]["Arithmetic Mean (Rang)"]
		uncorrectedStandardDeviation_1_Rang = stat_1["Dispersion"]["Uncorrected Standard Deviation (Rang)"]
		uncorrectedStandardDeviation_2_Rang = stat_2["Dispersion"]["Uncorrected Standard Deviation (Rang)"]
		
		results = {}
		
		def funcCovariance(sample_1, arithmeticMean_1, sample_2, arithmeticMean_2):
			result = 0
			for i in range(n):
				result += (sample_1[i]- arithmeticMean_1) * (sample_2[i] - arithmeticMean_2)
			result /= n
			return result
		results["Covariance (Rang)"] = rangedCovariance = funcCovariance(sampleRanged_1, arithmeticMean_1_Rang, sampleRanged_2, arithmeticMean_2_Rang)
		
		def funcSpearmansCorrelationCoefficient():
			result = rangedCovariance / (uncorrectedStandardDeviation_1_Rang * uncorrectedStandardDeviation_2_Rang)
			return result
		results["Spearman's Rang Correlation Coefficient"] = funcSpearmansCorrelationCoefficient()
		
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
		(index, sample, stat) = self.__params
		plt.ioff()	
		
		sns.set(color_codes=True)

		fig = plt.figure(figsize=(6, 6))

		number = stat["Binning"]["Number"]
		x_min = stat["Location"]["Min"] - (stat["Binning"]["Intervals"][1] - stat["Binning"]["Intervals"][0])/2
		x_max = stat["Location"]["Max"] + (stat["Binning"]["Intervals"][number] - stat["Binning"]["Intervals"][number-1])/2

		ax1 = plt.subplot2grid((15, 8), (0, 0), colspan=9)
		axBoxPlot = sns.boxplot(sample, ax=ax1, showfliers=False)

		ax2 = plt.subplot2grid((15, 8), (1, 0), colspan=9, rowspan=8)
		axDistPlot = sns.distplot(sample, ax=ax2, bins=stat["Binning"]["Intervals"], kde=False, norm_hist=False)

		ax3 = plt.subplot2grid((15, 8), (9, 0), colspan=9, rowspan=6)
		range = [x_min-(x_max-x_min)/100, x_max+(x_max-x_min)/100]
		Cde = plt.hist(sample, linewidth=2.0, range=range, histtype='step', cumulative=True, normed=1, bins=1000)

		x_limitsBoxPlot = axBoxPlot.set_xlim([x_min, x_max])
		x_limitsDistPlot = axDistPlot.set_xlim([x_min, x_max])
		CdeXlim = plt.xlim(x_min, x_max)
		CdeYlim = plt.ylim(-0.01, 1.01)

		axBoxXTickLabels = axBoxPlot.set_xticklabels([])
		axDistXTickLabels = axDistPlot.set_xticklabels([])
		axBoxYTicks = axBoxPlot.set_yticks([])
		axCdeYTicks = plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
	
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
		self.__correlations = {}
		
	
	@classmethod
	def create(cls, G):
		result = cls(G, cls.__TOKEN)
		
		for parameter in [ 
			(centrality.DegreeCentrality, 			(G, )),
			#(centrality.CoreDecomposition, 			(G, )),
			#(centrality.LocalClusteringCoefficient, (G, )),
			#(centrality.PageRank, 					(G, )),
			#(centrality.KPathCentrality,			(G, )),
			#(centrality.KatzCentrality,				(G, )),
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
		
		templateMeasure = readfile("measure.html")
		
		centralities = ""
		for key in self.__measures:
			measure = self.__measures[key]
			image = measure["image"]
			stat = measure["stat"]
			centralities += self.__formatMeasureTemplate(templateMeasure, key, image, stat)
			
			#centralities = centralities + "<div class=\"Plot\" title=\"" + key + "\" data-image=\"data:image/svg+xml;utf8," + image + "\" />"
			
		templateProfile = readfile("profile.html")
		result = templateProfile.format(**locals())
		display_html(HTML(result))
		
		if self.__verbose:
			print("\ntotal time: {:.2F} s".format(timerAll.elapsed))
	
		self.__pageCount = self.__pageCount + 1
	
	
	def __formatMeasureTemplate(self, template, key, image, stat):
		result = template.format(**locals())
		return result
	
	
	def __addMeasure(self, args):
		(measureClass, parameters) = args
		measureName = measureClass.__name__
		measure = {}
		measure["index"] = len(self.__measures)
		measure["class"] = measureClass
		measure["parameters"] = parameters
		measure["data"] = {}
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
			measure["data"]["sample"] = sample = instance.scores()
			measure["data"]["sorted"] = sampleSorted = sorted(measure["data"]["sample"])
			measure["data"]["ranged"] = sampleRanged = ranged(measure["data"]["sample"])
			tasks.put(Stat_Task(name, (sample, sampleSorted, sampleRanged)))
			numberOfTasks += 1
			measure["time"] = elapsed
		
		while(numberOfTasks):
			(type, name, data) = results.get()
			if (type == "Plot"):
				(index, image) = data
				self.__measures[name]["image"] = image
			elif (type == "Stat"):
				self.__measures[name]["stat"] = data
				tasks.put(Plot_Task(name, (0, self.__measures[name]["data"]["sample"], data)))
				numberOfTasks += 1
				
				for key in self.__correlations:
					print(name + " <-> " + key)
				self.__correlations[name] = {}
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