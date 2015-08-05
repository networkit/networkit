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
import collections
from array import array

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def readfile(postfix):
	""" private helper function: profiling-meta-file to string (all whitespace characters are replaced by ' ') """
	with open(__file__[:__file__.rfind(".py")] + "." + postfix, "r") as file:
		return " ".join(file.read().split())

		
try:
	__IPYTHON__
	
	def _initHeader(tag, type, data):
		""" private helper function for notebook hack: create content of extended header """
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

		
	def _initOverlay(name, data):
		""" private helper function for notebook hack: create content of overlay """
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
				""" + _initHeader("script", "javascript", readfile("js"))  + """
				""" + _initHeader("style",  "css",        readfile("css")) + """
				""" + _initOverlay("Overlay", readfile("overlay.html")) + """
			-->
			</script>
		""")
	)
except Exception as e:
	print(str(e))
	

def ranged(sample):
	""" TODO: """
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

	
class Theme:
	""" TODO: """
	def __init__(self):
		self.set()
    
	
	@classmethod
	def RGBA2RGB(cls, color, alpha, background):
		result = (
			color[0] * alpha + background[0] * (1-alpha),
			color[1] * alpha + background[1] * (1-alpha),
			color[2] * alpha + background[2] * (1-alpha),
			1
		)
		return result		
	
	
	def set(self, style="light", color=(0, 0, 1)):
		optionsStyle = ["light"]
		if style not in optionsStyle:
			raise ValueError("possible style options: " + str(optionsStyle))
		if len(color) != 3:
			raise ValueError("(r,g,b) tuple required")
			
		if style == "light":
			self.__color = color
			self.__defaultColor = (0, 0, 0)
			self.__defaultWidth = 1
			self.__backgroundColor = (1, 1, 1)
			self.__plotColor = Theme.RGBA2RGB(color, 0.6, self.__backgroundColor)
			self.__plotWidth = 3
			self.__faceColor = (color[0], color[1], color[2], 0.2)
			self.__faceColorGray = "lightgray"
			self.__edgeColor = (color[0], color[1], color[2], 0.6)
			self.__edgeColorGray = (0, 0, 0)
			self.__edgeWidth = 2
			self.__gridColor = "lightgray"
			self.__fontColor = (0, 0, 0)
            
		self.__fontSize = 10
		self.__style = style
    
    
	def get(self):
		return (self.__style, self.__color)
		
		
	def getDefaultColor(self):
		return self.__defaultColor
    
    
	def getDefaultWidth(self):
		return self.__defaultWidth
    
    
	def getPlotColor(self):
		return self.__plotColor
    
    
	def getPlotWidth(self):
		return self.__plotWidth
      
    
	def getFaceColor(self):
		return self.__faceColor
    
	
	def getFaceColorGray(self):
		return self.__faceColorGray
    
    
	def getEdgeColor(self):
		return self.__edgeColor
    
	
	def getEdgeColorGray(self):
		return self.__edgeColorGray
    
    
	def getEdgeWidth(self):
		return self.__edgeWidth
    
    
	def getBackgroundColor(self):
		return self.__backgroundColor
    
    
	def getGridColor(self):
		return self.__gridColor
    
    
	def getFontSize(self):
		return self.__fontSize
    
    
	def getFontColor(self):
		return self.__fontColor
    
	
class Worker(multiprocessing.Process):
	""" TODO: """
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
			data = "Error"
			try:
				data = task.run()
			except Exception as e:
				print("Error: " + task.getType() + " - " + task.getName())
				print(str(e))
			result = (task.getType(), task.getName(), data)
			self.__tasks.task_done()
			self.__results.put(result)


class ThreadPool():
	""" TODO: """
	def __init__(self, numberOfWorkers):
		self.__numberOfWorkers = numberOfWorkers
		self.__numberOfTasks = 0
		self.__tasks = multiprocessing.JoinableQueue()
		self.__results = multiprocessing.Queue()
		self.__workers = [Worker(self.__tasks, self.__results) for i in range(self.__numberOfWorkers)]
		for w in self.__workers:
			w.deamon = True
			w.start()
	
	
	def numberOfTasks(self):
		return self.__numberOfTasks
	
	
	def put(self, task):
		self.__tasks.put(task)
		self.__numberOfTasks += 1
		
	
	def get(self):
		result = self.__results.get()
		self.__numberOfTasks -= 1
		return result;

	
	def join(self):
		for i in range(self.__numberOfWorkers):
			self.__tasks.put(None)
		self.__tasks.join()
		
		
class Stat_Task:
	""" TODO: """
	def __init__(self, name, params):
		self.__name = name
		self.__params = params

	def getName(self):
		return self.__name

	def getType(self):
		return "Stat"

	def run(self):
		(sample, sampleSorted, sampleRanged, calculatePie) = self.__params
		n = len(sample)

		results = {}
		results["Properties"] = {}
		results["Location"] =  {}
		results["Dispersion"] =  {}
		results["Shape"] =  {}
		results["Binning"] = {}
		results["Distribution"] = {}

		results["Properties"]["Size"] = n

		def funcMin():
			result = sampleSorted[0]
			return result
		results["Location"]["Min"] = min = funcMin()

		def funcMax():
			result = sampleSorted[n-1]
			return result
		results["Location"]["Max"] = max = funcMax()

		def funcBesselsCorrection():
			try:
				result = n / (n-1)
			except:
				result = float("nan")
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
		results["Location"]["Quadratic Mean"] = quadraticMean = hoelderMean(sample, 2)
		results["Location"]["Cubic Mean"] = cubicMean = hoelderMean(sample, 3)
		if min > 0:
			results["Location"]["Harmonic Mean"] = harmonicMean = hoelderMean(sample, -1)
		else:
			results["Location"]["Harmonic Mean"] = harmonicMean = float("nan")

		def funcArithmeticMeanRang():
			result = (n + 1) / 2
			return result
		results["Location"]["Arithmetic Mean (Rang)"] = arithmeticMean_Rang = funcArithmeticMeanRang()


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
			result = (min + max)/ 2
			return result
		results["Location"]["Mid-Range"] = midRange = funcMidRange()

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
			return int(result)
		results["Binning"]["Number Histogram"] = k_Bins_Histogram = funcNumberOfBins()
		k_Bins_CDF = 256

		def funcIntervals(numberOfBins):
			result = []
			w = sampleRange / numberOfBins
			result.append(min)
			for i in range(1, numberOfBins):
				result.append(min + w * i)
			result.append(max)
			return result
		results["Binning"]["Intervals Histogram"] = intervalsHistogram = funcIntervals(k_Bins_Histogram)
		intervalsCDF = funcIntervals(k_Bins_CDF)
				
		def funcBinAbsoluteFrequencies(numberOfBins, intervals, comulative):
			result = []
			index = 0
			for i in range(numberOfBins):
				result.append(0)
			for i in range(n):
				value = sampleSorted[i]
				while intervals[index + 1] < value:
					index += 1
					if comulative and index > 0:
						result[index] = result[index-1]
				result[index] += 1
			return result
		results["Binning"]["Absolute Frequencies Histogram"] = absoluteFrequenciesHistogram = funcBinAbsoluteFrequencies(k_Bins_Histogram, intervalsHistogram, False)
		absoluteFrequenciesCDF = funcBinAbsoluteFrequencies(k_Bins_CDF, intervalsCDF, True)
		
		def funcJoinEmptyBins(k_Bin, intervals, frequencies, commulative):
			result = k_Bin
			value = 0
			for i in range(k_Bin):
				if frequencies[k_Bin-i-1] == value:
					del frequencies[k_Bin-i-1]
					del intervals[k_Bin-i]
					result -= 1
				if commulative:
					value = frequencies[k_Bin-i-1]
			return result
		results["Binning"]["Number CDF"] = k_Bins_CDF = funcJoinEmptyBins(k_Bins_CDF, intervalsCDF, absoluteFrequenciesCDF, True)
		results["Binning"]["Absolute Frequencies CDF"] = absoluteFrequenciesCDF
		results["Binning"]["Intervals CDF"] = intervalsCDF
		
		def funcBinRelativeFrequencies(absoluteFrequencies):
			result = []
			for H in absoluteFrequencies:
				result.append(H / n)
			return result
		results["Binning"]["Relative Frequencies Histogram"] = relativeFrequenciesHistogram = funcBinRelativeFrequencies(absoluteFrequenciesHistogram)
		results["Binning"]["Relative Frequencies CDF"] = relativeFrequenciesCDF = funcBinRelativeFrequencies(absoluteFrequenciesCDF)

		def funcMode():
			index = 0
			max = 0
			for i in range(len(absoluteFrequenciesHistogram)):
				if absoluteFrequenciesHistogram[i] > max:
					max = absoluteFrequenciesHistogram[i]
					index = i
			result = ((intervalsHistogram[index]+intervalsHistogram[index+1]) / 2, max)
			return result
		results["Binning"]["Mode"] = mode = funcMode()
		
		def funcLowerOutliers():
			lowerBound = Q1 - IQR * 3
			upperBound = Q1 - IQR * 1.5
			result_lower = min
			result_upper = min
			state = 0
			for i in range(n):
				value = sampleSorted[i]
				if value >= lowerBound and state == 0:
					result_lower = value 
					state = 1
				if value >= upperBound and state == 1:
					result_upper = value
					break
			return (result_lower, result_upper)
		results["Location"]["Outlier (Lower)"] = funcLowerOutliers()
		
		def funcUpperOutliers():
			lowerBound = Q3 + IQR * 1.5
			upperBound = Q3 + IQR * 3
			result_lower = max
			result_upper = max
			state = 0
			for i in range(n):
				value = sampleSorted[n-i-1]
				if value <= upperBound and state == 0:
					result_upper = value
					state = 1
				if value <= lowerBound and state == 1:
					result_lower = value
					break
			return (result_upper, result_lower)
		results["Location"]["Outlier (Upper)"] = funcUpperOutliers()
			
		if calculatePie:
			def funcPie():	
				n = len(sample)
				sum = 0
				for i in range(n):
					sum += sample[i]

				min = 0.015
				cutSize = 0
				cutValue = 0
				relativeFrequencies = [0]
				for i in range(n):
					value = sampleSorted[i]/sum
					if value < min:
						relativeFrequencies[0] += value
						cutSize += 1
						cutValue += value
					else:
						relativeFrequencies.append(value)
				return (relativeFrequencies, cutSize)
			results["Binning"]["Pie"] = funcPie()
		
			
		# Chi-Squared-Test <- Correct Binning
		# For Test-Case Purpose
		# n = 100
		# arithmeticMean = 51.05
		# s_n = 1.209
		# absoluteFrequencies = [5, 11, 35, 29, 13, 7]
		# intervals = [-10000, 49, 50, 51, 52, 53, 10000]
		# k_Bins = len(absoluteFrequencies)

		# def funcErf(x):
			# sign = 1 if x >= 0 else -1
			# x = abs(x)

			# a1 =  0.254829592
			# a2 = -0.284496736
			# a3 =  1.421413741
			# a4 = -1.453152027
			# a5 =  1.061405429
			# p  =  0.3275911

			# t = 1.0/(1.0 + p*x)
			# y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1) * t * math.exp(-x*x)
			# return sign*y

		# def funcDistributionNormal(x):
			# result = 1/2 * (1 + funcErf((x-arithmeticMean)/(math.sqrt(2) * s_n)))
			# return result

		# def funcDistributionExponential(x):
			# result = 1 - math.exp((-1/arithmeticMean) * x)
			# return result

		# def funcIncompleteGamma(s, x):
			# if x < 0.0:
				# return 0.0
			# sc = (1.0 / s)
			# sc *= math.pow(x, s)
			# sc *= math.exp(-x)
			# sum = 1.0
			# nom = 1.0
			# denom = 1.0
			# for i in range(1, 128):
				# nom *= x
				# denom *= s + i
				# sum += (nom / denom)
			# return sum * sc

		# def funcGamma(x):
			# result = (x / math.e) ** x
			# result *= math.sqrt(2 * math.pi / x)
			# result *= (1 + 1/(12 * x*x - 1/10)) ** x
			# return result

		# def funcPValue(criticalValue, degreesOfFreedom):
			# if criticalValue < 0.0 or degreesOfFreedom < 1:
				# return 0.0
			# k = degreesOfFreedom * 0.5
			# x = criticalValue * 0.5
			# if degreesOfFreedom == 2:
				# return math.exp(-x)
			# result = funcIncompleteGamma(k, x)
			# result /= funcGamma(k)
			# return 1-result

		# def funcChiSquaredTest(distribution, numberOfUsedEstimators):
			# z = 0;
			# pValue = 0
			# for i in range(k_Bins):
				# p_i = distribution(intervals[i+1]) - distribution(intervals[i])
				# hypotheticAbsoluteFrequency = n*p_i
				# if hypotheticAbsoluteFrequency == 0:
					# if absoluteFrequencies[i] == 0:
						# continue
					# z = float("Inf")
					# break
				# d = absoluteFrequencies[i] - hypotheticAbsoluteFrequency
				# z += d*d / hypotheticAbsoluteFrequency
				# print(self.getName(), hypotheticAbsoluteFrequency, absoluteFrequencies[i], z)
			# degreesOfFreedom = (k_Bins - 1) - numberOfUsedEstimators
			# pValue = funcPValue(z, degreesOfFreedom)
			# return (z, pValue)
		# results["Distribution"]["Chi-Square-Test (Normal)"] = funcChiSquaredTest(funcDistributionNormal, 2)
		# results["Distribution"]["Chi-Square-Test (Exponential)"] = funcChiSquaredTest(funcDistributionExponential, 1)

		return results


class Correlation_Task:
	""" TODO: """
	def __init__(self, name, params):
		self.__name = name
		self.__params = params

	def getName(self):
		return self.__name

	def getType(self):
		return "Correlation"

	def run(self):
		(nameB, sample_1, sampleRanged_1, stat_1, sample_2, sampleRanged_2, stat_2) = self.__params
		n = len(sample_1)
		assert (n == len(sample_2)), "sample sizes are not equal"

		results = {}

		def funcCovariance(sample_1, arithmeticMean_1, sample_2, arithmeticMean_2):
			result = 0
			for i in range(n):
				result += (sample_1[i]- arithmeticMean_1) * (sample_2[i] - arithmeticMean_2)
			result /= n
			return result
		results["Covariance"] = covariance = funcCovariance(
			sample_1,
			stat_1["Location"]["Arithmetic Mean"],
			sample_2,
			stat_2["Location"]["Arithmetic Mean"]
		)
		results["Covariance (Rang)"] = covarianceRanged = funcCovariance(
			sampleRanged_1,
			stat_1["Location"]["Arithmetic Mean (Rang)"],
			sampleRanged_2,
			stat_2["Location"]["Arithmetic Mean (Rang)"]
		)

		def funcPearsonsCorrelationCoefficient(covariance, uncorrectedStandardDeviation_1, uncorrectedStandardDeviation_2):
			result = covariance / (uncorrectedStandardDeviation_1 * uncorrectedStandardDeviation_2)
			return result
		results["Pearson's Correlation Coefficient"] = funcPearsonsCorrelationCoefficient(
			covariance,
			stat_1["Dispersion"]["Uncorrected Standard Deviation"],
			stat_2["Dispersion"]["Uncorrected Standard Deviation"]
		)
		results["Spearman's Rang Correlation Coefficient"] = funcPearsonsCorrelationCoefficient(
			covarianceRanged,
			stat_1["Dispersion"]["Uncorrected Standard Deviation (Rang)"],
			stat_2["Dispersion"]["Uncorrected Standard Deviation (Rang)"]
		)

		def funcFechnersCorrelationCoefficent(arithmeticMean_1, arithmeticMean_2):
			result = 0
			for i in range(n):
				result += math.copysign(1.0, (sample_1[i] - arithmeticMean_1) * (sample_2[i] - arithmeticMean_2))
			result /= n
			return result
		results["Fechner's Correlation Coefficient"] = funcFechnersCorrelationCoefficent(
			stat_1["Location"]["Arithmetic Mean"],
			stat_2["Location"]["Arithmetic Mean"]
		)

		return (nameB, results)


class PlotMeasure_Task:
	""" TODO: """
	def __init__(self, name, params):
		self.__name = name
		self.__params = params

	def getName(self):
		return self.__name

	def getType(self):
		return "PlotMeasure"

	def run(self):
		(index, stat, theme) = self.__params
		plt.ioff()

		
		def funcSpace(min, max):
			result = (max - min) * 0.04
			return result


		def funcTicks(min, max, numberOfTicks):
			result = []
			if numberOfTicks > 0:
				value = min    
				step = (max - min) / numberOfTicks
				for i in range(numberOfTicks):
					result.append(min + step*i)
					value += step
				result.append(max)
			return result


		def funcPlotEnd(fig, ax, theme, width, height, drawAxis=True):
			ax.patch.set_facecolor(theme.getBackgroundColor())
			if drawAxis:
				axisColor = theme.getGridColor()
			else:
				axisColor = theme.getBackgroundColor()
			ax.spines["bottom"].set_color(axisColor)
			ax.spines["top"].set_color(axisColor) 
			ax.spines["right"].set_color(axisColor)
			ax.spines["left"].set_color(axisColor)
			ax.tick_params(axis="x", colors=theme.getGridColor(), which="both", labelsize=theme.getFontSize())
			ax.tick_params(axis="y", colors=theme.getGridColor(), which="both", labelsize=theme.getFontSize())
			ax.xaxis.label.set_color(theme.getFontColor())
			ax.yaxis.label.set_color(theme.getFontColor())
			[x_ticklabel.set_color(theme.getFontColor()) for x_ticklabel in ax.get_xticklabels()]
			[y_ticklabel.set_color(theme.getFontColor()) for y_ticklabel in ax.get_yticklabels()]
			fig.set_size_inches(width, height)
			

		def funcPlotBox(ax, x_numberOfTicks, x_showTickLabels, showGrid):
			q1 = stat["Location"]["1st Quartile"]
			q3 = stat["Location"]["3rd Quartile"]
			median = stat["Location"]["Median"]
			outlier_lower = stat["Location"]["Outlier (Lower)"][0]
			outlier_upper = stat["Location"]["Outlier (Upper)"][0]
			whisker_lower = stat["Location"]["Outlier (Lower)"][1]
			whisker_upper = stat["Location"]["Outlier (Upper)"][1]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			space = funcSpace(x_min, x_max)
			ticks = funcTicks(x_min, x_max, x_numberOfTicks)
			ax.scatter(
				[x_min, x_max],
				[0.5, 0.5],
				color = theme.getEdgeColor(),
				s = 35
			)
			ax.scatter(
				[outlier_lower, outlier_upper],
				[0.5, 0.5],
				color = theme.getDefaultColor(),
				marker = 'x',
				s = 50
			)
			ax.plot(
				[median, median],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_lower, whisker_lower],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_lower, q1],
				[0.5, 0.5],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_upper, whisker_upper],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_upper, q3],
				[0.5, 0.5],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[median, median],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.add_patch(patches.Rectangle(
				(q1, 0.2),
				width = q3-q1,
				height = 0.6,
				facecolor = theme.getFaceColor(),
				linestyle = "solid",
				linewidth = theme.getEdgeWidth(),
				edgecolor = theme.getEdgeColor()
			))			
			ax.set_xlim([x_min-space, x_max+space])
			ax.set_ylim([0, 1])
			ax.set_xticks(ticks)
			ax.set_yticks([])
			if not x_showTickLabels:
				ax.set_xticklabels([])
			if showGrid:
				ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
			return ax
			

		def funcPlotHistogram(ax, x_numberOfTicks, y_numberOfTicks, x_showTickLabels, y_showTickLabels, showGrid):
			numberOfBins = stat["Binning"]["Number Histogram"]
			intervals = stat["Binning"]["Intervals Histogram"]
			absoluteFrequencies = stat["Binning"]["Absolute Frequencies Histogram"]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			y_min = 0
			y_max = stat["Binning"]["Mode"][1]
			x_space = funcSpace(x_min, x_max)
			y_space = funcSpace(y_min, y_max)
			x_ticks = funcTicks(x_min, x_max, x_numberOfTicks)
			for i in range(numberOfBins):
				ax.add_patch(patches.Rectangle(
					(intervals[i], 0),
					width = intervals[i+1]-intervals[i],
					height = absoluteFrequencies[i],
					facecolor = theme.getFaceColor(),
					linestyle = "solid",
					linewidth = theme.getEdgeWidth(),
					edgecolor = theme.getEdgeColor(),
				))
			ax.set_xlim([x_min-x_space, x_max+x_space])
			ax.set_ylim([0, y_max+y_space])
			ax.set_xticks(x_ticks)
			if not x_showTickLabels:
				ax.set_xticklabels([])
			if not y_showTickLabels:
				ax.set_yticklabels([])
			if showGrid:
				ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
			return ax


		def funcPlotCDE(ax, x_numberOfTicks, y_numberOfTicks, x_showTickLabels, y_showTickLabels, showGrid):
			numberOfBins = stat["Binning"]["Number CDF"]
			intervals = stat["Binning"]["Intervals CDF"]
			comulativeRelativeFrequencies = stat["Binning"]["Relative Frequencies CDF"]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			x_space = funcSpace(x_min, x_max)
			y_space = funcSpace(0, 1)
			x_ticks = funcTicks(x_min, x_max, x_numberOfTicks)
			y_ticks = funcTicks(0, 1, y_numberOfTicks)
			for i in range(numberOfBins):
				ax.plot(
					[intervals[i], intervals[i+1]],
					[comulativeRelativeFrequencies[i], comulativeRelativeFrequencies[i]],
					color = theme.getPlotColor(),
					linestyle = "-",
					linewidth = theme.getEdgeWidth()
				)
				ax.plot(
					[intervals[i], intervals[i]],
					[0 if i==0 else comulativeRelativeFrequencies[i-1], comulativeRelativeFrequencies[i]],
					color = theme.getDefaultColor(),
					linestyle = "dotted",
				linewidth = theme.getDefaultWidth()
				)
			ax.set_xlim([x_min-x_space, x_max+x_space])
			ax.set_ylim([0, 1+y_space])
			ax.set_xticks(x_ticks)
			ax.set_yticks(y_ticks)
			if not x_showTickLabels:
				ax.set_xticklabels([])
			if not y_showTickLabels:
				ax.set_yticklabels([])
			if showGrid:
				ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
			return ax
    
		def funcPlotPie(ax):
			numberOfTooSmallSubsets = stat["Binning"]["Pie"][1]
			relativeFrequencies = stat["Binning"]["Pie"][0]
			radius = 2
			accumulator = 0
			
			for i in range(len(relativeFrequencies)):
				value = relativeFrequencies[i]
				alpha = 360 * value
				label = "{:1.1f}%".format(value * 100)
				labelRadius = radius * 1.1
				if i == 0:
					label = str(numberOfTooSmallSubsets) + " Subsets\n" + label
				else:
					t = accumulator + alpha/2
					if (value < 0.022 and (
						t <= 30 and i%2 == 0 or
						t >= 150 and t <= 180 and i%2 == 1 or
						t >= 180 and t <= 210 and i%2 == 0 or
						t >= 330 and i%2 == 1
					)):
						labelRadius = radius * 1.19
				if accumulator + alpha/2 > 180:
					ha = "left"
				else:
					ha = "right"
				
				if i == 0:
					ax.add_patch(patches.Wedge(
						(math.cos(math.pi/180 * (90 + alpha/2)) * radius * 0.1, 
						math.sin(math.pi/180 * (90 + alpha/2)) * radius * 0.1),
						radius,
						90,
						90 + accumulator + alpha,
						facecolor = theme.getFaceColorGray(),
						edgecolor = theme.getDefaultColor()
					))
					labelRadius = radius * 1.18
				else:
					scale = 1-relativeFrequencies[1]/value
					faceColor = (
						theme.getPlotColor()[0] * scale,
						theme.getPlotColor()[1] * scale,
						theme.getPlotColor()[2] * scale,
						theme.getEdgeColor()[3]
					)
					ax.add_patch(patches.Wedge(
						(0, 0),
						radius,
						90 + accumulator,
						90 + accumulator + alpha,
						facecolor = faceColor,
						edgecolor = theme.getDefaultColor()
					))
				plt.text(
					math.cos(math.pi/180 * (90 + accumulator + alpha/2)) * labelRadius, 
					math.sin(math.pi/180 * (90 + accumulator + alpha/2)) * labelRadius,
					s = label,
					ha = ha,
					#family = 'sans-serif',
					size = theme.getFontSize()
				)
				accumulator += alpha
			
			ax.set_xlim([-3.2, 3.2])
			ax.set_ylim([-2.7, 2.7])
			ax.set_xticks([])
			ax.set_yticks([])
			return ax
	
	
		if index == 0:
			fig = plt.figure()
			
			ax1 = plt.subplot2grid((40, 8), (0, 0), colspan=8, rowspan=3)
			funcPlotBox(
				ax = ax1,
				x_numberOfTicks = 5,
				x_showTickLabels = False,
				showGrid = False
			)
			funcPlotEnd(
				fig = fig,
				ax = ax1,
				theme = theme,
				width = 4,
				height = 0.2
			)
			
			ax2 = plt.subplot2grid((40, 8), (3, 0), colspan=8, rowspan=20)
			funcPlotHistogram(
				ax = ax2,
				x_numberOfTicks = 5,
				y_numberOfTicks = 5,
				x_showTickLabels = False,
				y_showTickLabels = True,
				showGrid = True
			)
			funcPlotEnd(
				fig = fig,
				ax = ax2,
				theme = theme,
				width = 4,
				height = 3
			)

			ax3 = plt.subplot2grid((40, 8), (23, 0), colspan=8, rowspan=17)
			funcPlotCDE(
				ax = ax3,
				x_numberOfTicks = 5,
				y_numberOfTicks = 5,
				x_showTickLabels = True,
				y_showTickLabels = True,
				showGrid = True
			)
			funcPlotEnd(
				fig = fig,
				ax = ax3,
				theme = theme,
				width = 4,
				height = 3
			)
			fig.set_size_inches(6, 6)
		
		elif index == 1:
			fig = plt.figure()
			ax = fig.gca()
			
			funcPlotHistogram(
				ax = ax,
				x_numberOfTicks = 5,
				y_numberOfTicks = 5,
				x_showTickLabels = True,
				y_showTickLabels = True,
				showGrid = True
			)
			funcPlotEnd(
				fig = fig,
				ax = ax,
				theme = theme,
				width = 4,
				height = 2.5
			)
		
		elif index == 2:
			fig = plt.figure()
			ax = fig.gca()
			
			funcPlotPie(
				ax = ax,
			)
			funcPlotEnd(
				fig = fig,
				ax = ax,
				theme = theme,
				width = 6.4*1.5,
				height = 5.4*1.5,
				drawAxis = False
			)
	
		fig.tight_layout()
		imgdata = io.StringIO()
		fig.savefig(imgdata, format='svg')
		plt.close(fig)
		plaintext = imgdata.getvalue()
		plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
		encoded = quote(plaintext, safe='');
		return (index, encoded)


class PlotCorrelation_Task:
	""" TODO: """
	def __init__(self, name, params):
		self.__name = name
		self.__params = params

	def getName(self):
		return self.__name

	def getType(self):
		return "PlotCorrelation"

	def run(self):
		nameA = self.__name
		(nameB, sample_1, sample_2) = self.__params
		plt.ioff()

		fig = plt.figure(figsize=(4, 3.75))

		def hexbin(x, y, color, **kwargs):
			# cmap = sns.light_palette(color, as_cmap=True)
			ax = plt.hexbin(x, y, gridsize=32, bins="log", **kwargs)
			# ax = plt.hexbin(x, y, gridsize=32, bins="log", cmap=cmap, **kwargs)
			xLabel = plt.xlabel(nameA)
			yLabel = plt.ylabel(nameB)

		hexbin(sample_1, sample_2, "#000070")

		fig.tight_layout()
		imgdata = io.StringIO()
		fig.savefig(imgdata, format='svg')
		plt.close(fig)
		plaintext = imgdata.getvalue()
		plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
		encoded = quote(plaintext, safe='');
		return (nameB, encoded)


class Profile:
	""" TODO: """
	__TOKEN = object();
	__pageCount = 0
	__verbose = False
	__verboseLevel = 0
	__parallel = multiprocessing.cpu_count() * 2


	def __init__(self, G, token=object()):
		""" TODO: """
		if token is not self.__TOKEN:
			raise ValueError("call create(G) to create an instance")
		self.__G = G
		self.__properties = {}
		self.__measures = collections.OrderedDict()
		self.__correlations = {}


	@classmethod
	def create(cls, G, exclude=["KPathCentrality", "KatzCentrality"]):
		""" TODO: """
		result = cls(G, cls.__TOKEN)

		def funcScores(instance):
			return instance.scores()

		def funcSizes(instance):
			return sorted(instance.getPartition().subsetSizes())

		if G.isDirected():
			classConnectedComponents = properties.StronglyConnectedComponents
		else:
			classConnectedComponents = properties.ConnectedComponents

		for parameter in [
			("Node Centrality",	True,	funcScores,	centrality.DegreeCentrality, 			(G, )),
			("Node Centrality",	True,	funcScores,	centrality.CoreDecomposition, 			(G, )),
			("Node Centrality",	True,	funcScores,	centrality.LocalClusteringCoefficient,	(G, )),
			("Node Centrality",	True,	funcScores,	centrality.PageRank, 					(G, )),
			("Node Centrality",	True,	funcScores,	centrality.KPathCentrality,				(G, )),
			("Node Centrality",	True,	funcScores,	centrality.KatzCentrality,				(G, )),
			("Node Centrality",	True,	funcScores,	centrality.ApproxBetweenness2,			(G, max(42, G.numberOfNodes() / 10000), False)),
		#	("Partition",		False,	funcSizes,	community.LPDegreeOrdered, 				(G, )),
		#	("Partition",		False,	funcSizes,	community.PLP, 							(G, )),
			("Partition",		False,	funcSizes,	community.PLM, 							(G, )),
			("Partition",		False,	funcSizes,	classConnectedComponents,				(G, ))
		]: result.__addMeasure(parameter, exclude)

		result.__loadProperties()
		result.__loadMeasures()
		return result;


	@classmethod
	def setVerbose(cls, verbose=False, level=0):
		""" TODO: """
		cls.__verbose = verbose
		cls.__verboseLevel = level


	@classmethod
	def getVerbose(cls):
		""" TODO: """
		return (cls.__verbose, cls.__verboseLevel)


	@classmethod
	def setParallel(cls, parallel):
		""" TODO: """
		if (parallel < 1):
			raise ValueError("parallel < 1");
		cls.__parallel = parallel


	@classmethod
	def getParallel(cls):
		""" TODO: """
		return cls.__parallel


	def getStat(self, measure):
		""" TODO: """
		return self.__measures[measure]["stat"]


	def getCategory(self, measure):
		""" TODO: """
		return self.__measures[measure]["category"]


	def getElapsedTime(self, measure):
		""" TODO: """
		return self.__measures[measure]["time"]


	def show(self, theme=Theme()):
		""" TODO: """
		try:
			__IPYTHON__
		except:
			raise RuntimeError("this function cannot be used outside ipython notebook")
			
		if self.__verbose:
			timerAll = stopwatch.Timer()
		
		pool = ThreadPool(self.__parallel)
		for name in self.__measures:
			category = self.__measures[name]["category"]
			pool.put(
				PlotMeasure_Task(name, (
					0,
					self.__measures[name]["stat"],
					theme
				))
			)
			pool.put(
				PlotMeasure_Task(name, (
					1,
					self.__measures[name]["stat"],
					theme
				))
			)				
			if category == "Partition":
				pool.put(
					PlotMeasure_Task(name, (
						2,
						self.__measures[name]["stat"],
						theme
					))
				)
		while pool.numberOfTasks() > 0:
			(type, name, data) = pool.get()
			category = self.__measures[name]["category"]

			if type == "PlotMeasure":
				(index, image) = data
				self.__measures[name]["image"][index] = image	
		pool.join()

		templateMeasure = readfile("measure.html")

		results = {}
		for category in self.__correlations:
			results[category] = {}
			results[category]["Correlations"] = {}
			results[category]["Correlations"]["HeatMaps"] = ""
			results[category]["Correlations"]["ScatterPlots"] = ""
			results[category]["Measures"] = ""
			results[category]["Overview"] = ""

			def funcHeatMap(category, correlationName):
				result = "<div class=\"SubCategory HeatTable\" data-title=\"" + correlationName + "\">"
				keyBList = []
				for keyA in self.__measures:
					if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
						keyBList.append(keyA)
						for keyB in keyBList:
							try:
								value = self.__correlations[category][keyA][keyB]
							except:
								value = self.__correlations[category][keyB][keyA]
							result += "<div class=\"HeatCell\" title=\"" + keyB + " - " + keyA + "\" data-image=\"data:image/svg+xml;utf8," + value["image"] + "\" data-heat=\"{:+.3F}\"></div>".format(value["stat"][correlationName])
						result += "<div class=\"HeatCellName\">" + keyB + "</div><br>"
				result += "</div>"
				return result
			results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Pearson's Correlation Coefficient")
			results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Spearman's Rang Correlation Coefficient")
			results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Fechner's Correlation Coefficient")

			def funcScatterPlot(category):
				result = ""
				keyBList = []
				for keyA in self.__measures:
					if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
						keyBList.append(keyA)
						for keyB in keyBList:
							if keyA != keyB:
								try:
									value = self.__correlations[category][keyA][keyB]
								except:
									value = self.__correlations[category][keyB][keyA]
								result += "<div class=\"Thumbnail_ScatterPlot\" data-title=\"" + keyB + "\" data-title-second=\"" + keyA + "\"><img src=\"data:image/svg+xml;utf8," + value["image"] + "\" /></div>"
				return result
			results[category]["Correlations"]["ScatterPlots"] += funcScatterPlot(category)

		for key in self.__measures:
			measure = self.__measures[key]
			category = measure["category"]
			image = measure["image"]
			stat = measure["stat"]
			results[category]["Measures"] += self.__formatMeasureTemplate(
				templateMeasure,
				key,
				image,
				stat
			)
			try:
				results[category]["Measures"] += "<div class=\"PartitionPie\"><img src=\"data:image/svg+xml;utf8," + image[2] + "\" /></div>"
			except:
				pass
			results[category]["Overview"] += "<div class=\"Thumbnail_Overview\" data-title=\"" + key + "\"><a href=\"#NetworKit_Page_" + str(self.__pageCount) + "_" + key + "\"><img src=\"data:image/svg+xml;utf8," + image[1] + "\" /></a></div>"

		templateProfile = readfile("profile.html")
		result = self.__formatProfileTemplate(
			templateProfile,
			results
		)
		display_html(HTML(result))
		self.__pageCount = self.__pageCount + 1

		if self.__verbose:
			print("\ntotal time: {:.2F} s".format(timerAll.elapsed))


	def __formatMeasureTemplate(self, template, key, image, stat):
		""" TODO: """
		pageIndex = self.__pageCount
		result = template.format(**locals())
		return result


	def __formatProfileTemplate(self, template, results):
		""" TODO: """
		pageIndex = self.__pageCount
		properties = self.__properties
		result = template.format(**locals())
		return result


	def __addMeasure(self, args, exclude):
		""" TODO: """
		(measureCategory, correlate, getter, measureClass, parameters) = args
		measureName = measureClass.__name__
		if measureName not in exclude:
			measure = {}
			measure["category"] = measureCategory
			measure["correlate"] = correlate
			measure["getter"] = getter
			measure["class"] = measureClass
			measure["parameters"] = parameters
			measure["data"] = {}
			measure["image"] = {}
			self.__measures[measureName] = measure
		try:
			self.__correlations[measureCategory]
		except:
			self.__correlations[measureCategory] = {}


	def __loadProperties(self):
		""" TODO: """
		self.__properties["Nodes"] = self.__G.numberOfNodes()
		self.__properties["Edges"] = self.__G.numberOfEdges()
		self.__properties["Directed"] = self.__G.isDirected()
		self.__properties["Weighted"] = self.__G.isWeighted()
		self.__properties["Density"] = properties.density(self.__G)
		self.__properties["Diameter Range"] = properties.Diameter.estimatedDiameterRange(self.__G, error=0.1)


	def __loadMeasures(self):
		""" TODO: """
		def funcPrint(str):
			if self.__verbose:
				if self.__verboseLevel >= 1:
					print(str, flush=True)
				else:
					print(".", end="", flush=True)

		pool = ThreadPool(self.__parallel)

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
				print("{:.2F} s".format(elapsed), flush=True)
			measure["data"]["sample"] = measure["getter"](instance)
			measure["data"]["sorted"] = sorted(measure["data"]["sample"])
			measure["data"]["ranged"] = ranged(measure["data"]["sample"])
			measure["time"] = elapsed

		if self.__verbose:
			print("")

		for name in self.__measures:
			if len(self.__measures[name]["data"]["sample"]) <= 1:
				del self.__measures[name]
			else:
				category = self.__measures[name]["category"]
				pool.put(
					Stat_Task(name, (
						self.__measures[name]["data"]["sample"],
						self.__measures[name]["data"]["sorted"],
						self.__measures[name]["data"]["ranged"],
						category == "Partition"
					))
				)
				
		while pool.numberOfTasks() > 0:
			(type, name, data) = pool.get()
			category = self.__measures[name]["category"]

			if type == "Stat":
				self.__measures[name]["stat"] = data
				funcPrint("Stat: " + name)					
				if self.__measures[name]["correlate"]:
					for key in self.__correlations[category]:
						self.__correlations[category][key][name] = {}
						self.__correlations[category][key][name]["stat"] = {}
						pool.put(
							Correlation_Task(key, (
								name,
								self.__measures[key]["data"]["sample"],
								self.__measures[key]["data"]["ranged"],
								self.__measures[key]["stat"],
								self.__measures[name]["data"]["sample"],
								self.__measures[name]["data"]["ranged"],
								self.__measures[name]["stat"]
							))
						)
						
						pool.put(
							PlotCorrelation_Task(key, (
								name,
								self.__measures[key]["data"]["sample"],
								self.__measures[name]["data"]["sample"]
							))
						)
						
					self.__correlations[category][name] = {}
					self.__correlations[category][name][name] = {}
					self.__correlations[category][name][name]["stat"] = {
						"Spearman's Rang Correlation Coefficient": 1,
						"Pearson's Correlation Coefficient": 1,
						"Fechner's Correlation Coefficient": 1
					}
					self.__correlations[category][name][name]["image"] = ""

			elif type == "Correlation":
				(nameB, correlation) = data
				funcPrint("Correlation: " + name + " <-> " + nameB)
				self.__correlations[category][name][nameB]["stat"] = correlation

			elif type == "PlotCorrelation":
				(nameB, image) = data
				funcPrint("Plot (Correlation): " + name)
				self.__correlations[category][name][nameB]["image"] = image

		pool.join()

		if self.__verbose:
			if self.__verboseLevel < 1:
				print("")
			print("\ntotal time (measures + stats + correlations + plots): {:.2F} s".format(timerAll.elapsed))