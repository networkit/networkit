#
# file: stat.py
# author: Mark Erb
#

from . import job

import math
import matplotlib.pyplot as plt

from _NetworKit import sort2
from _NetworKit import ranked


def sorted(sample):
	"""	returns a sorted list of given numbers """
	return sort2(sample)

class Stat(job.Job):
	""" statistical computation object """
	
	def __init__(self, name, params):
		""" constructor: see PlotJob and .run() """
		job.Job.__init__(
			self,
			"Stat",
			name
		)
		self.__params = params

	def run(self):
		""" computation """
		(sample, sampleSorted, sampleRanked, calculatePie) = self.__params
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

		def funcArithmeticMeanRank():
			result = (n + 1) / 2
			return result
		results["Location"]["Arithmetic Mean (Rank)"] = arithmeticMean_Rank = funcArithmeticMeanRank()

		def funcUncorrectedVariance(sample, arithmeticMean):
			result = 0
			for i in range(n):
				result += (sample[i] - arithmeticMean) ** 2
			result /= n
			return result
		results["Dispersion"]["Uncorrected Variance"] = variance_uncorrected = funcUncorrectedVariance(sample, arithmeticMean)
		results["Dispersion"]["Uncorrected Variance (Rank)"] = variance_Rank_uncorrected = funcUncorrectedVariance(sampleRanked, arithmeticMean_Rank)

		def funcVariance(variance_uncorrected):
			result = variance_uncorrected * besselsCorrection
			return result
		results["Dispersion"]["Variance"] = variance = funcVariance(variance_uncorrected)
		results["Dispersion"]["Variance (Rank)"] = variance_Rank = funcVariance(variance_Rank_uncorrected)

		def funcStandardDeviation(variance):
			result = math.sqrt(variance)
			return result
		results["Dispersion"]["Standard Deviation"] = s_n = funcStandardDeviation(variance)
		results["Dispersion"]["Standard Deviation (Rank)"] = s_n_Rank = funcStandardDeviation(variance_Rank)
		results["Dispersion"]["Uncorrected Standard Deviation"] = s_n_uncorrected = funcStandardDeviation(variance_uncorrected)
		results["Dispersion"]["Uncorrected Standard Deviation (Rank)"] = s_n_Rank_uncorrected = funcStandardDeviation(variance_Rank_uncorrected)

		def funcCoefficientOfVariation(s_n, arithmeticMean):
			result = float("nan")
			if arithmeticMean != 0:
				result = s_n / arithmeticMean
			return result
		results["Dispersion"]["Coefficient Of Variation"] = c_v = funcCoefficientOfVariation(s_n, arithmeticMean)
		results["Dispersion"]["Coefficient Of Variation (Rank)"] = c_v_Rank = funcCoefficientOfVariation(s_n_Rank, arithmeticMean_Rank)
		results["Dispersion"]["Uncorrected Coefficient Of Variation"] = c_v = funcCoefficientOfVariation(s_n_uncorrected, arithmeticMean)
		results["Dispersion"]["Uncorrected Coefficient Of Variation (Rank)"] = c_v_Rank = funcCoefficientOfVariation(s_n_Rank_uncorrected, arithmeticMean_Rank)
		
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
			result = float("nan")
			if s_n != 0:
				result = 3 * (arithmeticMean - median) / s_n
			return result
		results["Shape"]["Skewness YP"] = skewness_yp = funcSkewnessYP()

		def funcMomentum(p):
			result = float("nan")
			if s_n != 0:
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

		def funcNumberOfBins(commulative):
			result = 1
			if min < max:
				if commulative:
					result = 256
				else:
					result = math.sqrt(n)
					if (result < 5):
						result = 5
					elif (result > 20):
						result = 20
			return int(result)
		results["Binning"]["Number Histogram"] = k_Bins_Histogram = funcNumberOfBins(False)
		k_Bins_CDF = funcNumberOfBins(True)

		def funcIntervals(numberOfBins):
			result = []
			w = sampleRange / numberOfBins
			result.append(min)
			for i in range(1, numberOfBins):
				result.append(min + w * i)
			result.append(max if min < max else max+10e-12)
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


		# The following code within this class is experimental and under develop
		#
		# Chi-Squared-Test <- Correct Binning
		# For Test-Case Purpose uncomment the following lines 
		# n = 100
		# arithmeticMean = 51.05
		# s_n = 1.209
		# absoluteFrequencies = [5, 11, 35, 29, 13, 7]
		# intervals = [-10000, 49, 50, 51, 52, 53, 10000]
		# k_Bins = len(absoluteFrequencies)

		def funcErf(x):
			sign = 1 if x >= 0 else -1
			x = abs(x)

			a1 =  0.254829592
			a2 = -0.284496736
			a3 =  1.421413741
			a4 = -1.453152027
			a5 =  1.061405429
			p  =  0.3275911

			t = 1.0/(1.0 + p*x)
			y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1) * t * math.exp(-x*x)
			return sign*y

		def funcDistributionNormal(x):
			result = 1/2 * (1 + funcErf((x-arithmeticMean)/(math.sqrt(2) * s_n)))
			return result

		def funcDistributionExponential(x):
			result = 1 - math.exp((-1/arithmeticMean) * x)
			return result
			
		def funcDistributionExponentialInverse(x):
			result = math.ln(1/(1-a))*arithmeticMean
			return result
			
		def funcIncompleteGamma(s, x):
			if x < 0.0:
				return 0.0
			sc = (1.0 / s)
			sc *= math.pow(x, s)
			sc *= math.exp(-x)
			sum = 1.0
			nom = 1.0
			denom = 1.0
			for i in range(1, 128):
				nom *= x
				denom *= s + i
				sum += (nom / denom)
			return sum * sc

		def funcGamma(x):
			result = (x / math.e) ** x
			result *= math.sqrt(2 * math.pi / x)
			result *= (1 + 1/(12 * x*x - 1/10)) ** x
			return result

		def funcPValue(criticalValue, degreesOfFreedom):
			if criticalValue < 0.0 or degreesOfFreedom < 1:
				return 0.0
			k = degreesOfFreedom * 0.5
			x = criticalValue * 0.5
			if degreesOfFreedom == 2:
				return math.exp(-x)
			result = funcIncompleteGamma(k, x)
			result /= funcGamma(k)
			return 1-result

		def funcNumberOfBinsChiSquaredTest():
			result = 1 + ln(n)/ln(2)
			if result > 128:
				result = 128
			return int(result)
			
		# TODO:
		# def funcIntervalsChiSquaredTest(inverseFunction, min, max, k_Bins):
			# bin_size = n / k
			# result = []
			# return result
			
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


class Correlation(job.Job):
	""" correlation computation object """
	
	def __init__(self, name, params):
		""" constructor: see PlotJob and .run() """
		job.Job.__init__(
			self,
			"Correlation",
			name
		)
		self.__params = params

		
	def run(self):
		""" computation """
		(nameB, sample_1, sampleRanked_1, stat_1, sample_2, sampleRanked_2, stat_2) = self.__params
		n = len(sample_1)
		assert (n == len(sample_2)), "sample sizes are not equal"

		results = {}
		results["Value"] = {}

		def funcCovariance(sample_1, arithmeticMean_1, sample_2, arithmeticMean_2):
			result = 0
			for i in range(n):
				result += (sample_1[i]- arithmeticMean_1) * (sample_2[i] - arithmeticMean_2)
			result /= n
			return result
		results["Value"]["Covariance"] = covariance = funcCovariance(
			sample_1,
			stat_1["Location"]["Arithmetic Mean"],
			sample_2,
			stat_2["Location"]["Arithmetic Mean"]
		)
		results["Value"]["Covariance (Rank)"] = covarianceRanked = funcCovariance(
			sampleRanked_1,
			stat_1["Location"]["Arithmetic Mean (Rank)"],
			sampleRanked_2,
			stat_2["Location"]["Arithmetic Mean (Rank)"]
		)

		def funcPearsonsCorrelationCoefficient(covariance, uncorrectedStandardDeviation_1, uncorrectedStandardDeviation_2):
			result = float("nan")
			if uncorrectedStandardDeviation_1 * uncorrectedStandardDeviation_2 != 0:
				result = covariance / (uncorrectedStandardDeviation_1 * uncorrectedStandardDeviation_2)
			return result
		results["Value"]["Pearson's Correlation Coefficient"] = funcPearsonsCorrelationCoefficient(
			covariance,
			stat_1["Dispersion"]["Uncorrected Standard Deviation"],
			stat_2["Dispersion"]["Uncorrected Standard Deviation"]
		)
		results["Value"]["Spearman's Rank Correlation Coefficient"] = funcPearsonsCorrelationCoefficient(
			covarianceRanked,
			stat_1["Dispersion"]["Uncorrected Standard Deviation (Rank)"],
			stat_2["Dispersion"]["Uncorrected Standard Deviation (Rank)"]
		)

		def funcFechnersCorrelationCoefficent(arithmeticMean_1, arithmeticMean_2):
			result = 0
			for i in range(n):
				result += math.copysign(1.0, (sample_1[i] - arithmeticMean_1) * (sample_2[i] - arithmeticMean_2))
			result /= n
			return result
		results["Value"]["Fechner's Correlation Coefficient"] = funcFechnersCorrelationCoefficent(
			stat_1["Location"]["Arithmetic Mean"],
			stat_2["Location"]["Arithmetic Mean"]
		)
		
		def funcHexBinning(sample_1, sample_2):
			""" binning for scatter plots """
			result = {}
			n = 32
			
			fig = plt.figure()
			extent = [
				stat_1["Location"]["Min"],
				stat_1["Location"]["Max"],
				stat_2["Location"]["Min"],
				stat_2["Location"]["Max"]
			]
			image = plt.hexbin(
				sample_1,
				sample_2,
				gridsize = n,
				extent = extent
			)
			
			result["Grid Size"] = n
			result["Absolute Frequencies"] = frequencies = image.get_array()
			max = 0
			for value in frequencies:
				if max < value:
					max = value
			result["Max Frequency"] = max
			result["Offsets"] = image.get_offsets()
			result["Paths"] = image.get_paths()[0]
			plt.close(fig)
			return result
		results["Binning"] = funcHexBinning(
			sample_1,
			sample_2
		) 

		return (nameB, results)
