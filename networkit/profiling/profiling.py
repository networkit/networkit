#
# file: profiling.py
# author: Mark Erb
#

from networkit import *
import networkit as kit

from . import multiprocessing
from . import stat
from . import plot

from IPython.core.display import *
import collections
import math


def readfile(filename, removeWS=False):
	""" private helper function for file-loading """
	with open(__file__[:__file__.rfind("/")+1] + filename, "r") as file:
		result = file.read()
		if removeWS:
			result = " ".join(result.split())
		return result


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
				""" + _initHeader("script", "javascript", readfile("html/profiling.js", True))  + """
				""" + _initHeader("style",  "css",        readfile("html/profiling.css", True)) + """
				""" + _initOverlay("Overlay", readfile("html/overlay.html", True)) + """
			-->
			</script>
		""")
	)
except Exception as e:
	print(str(e))


class Profile:
	""" TODO: """
	__TOKEN = object();
	__pageCount = 0
	__verbose = False
	__verboseLevel = 0
	__verboseFilename = ""
	__parallel = multiprocessing.numberOfProcessors()


	def __init__(self, G, token=object()):
		""" TODO: """
		if token is not self.__TOKEN:
			raise ValueError("call create(G) to create an instance")
		self.__G = G
		self.__properties = {}
		self.__measures = collections.OrderedDict()
		self.__correlations = {}


	@classmethod
	def create(cls, G, exclude=["KPathCentrality"]):
		""" TODO: """
		result = cls(G, cls.__TOKEN)
		kit.setNumberOfThreads(result.__parallel)

		def funcScores(instance):
			return instance.scores()

		def funcSizes(instance):
			return sorted(instance.getPartition().subsetSizes())



		if G.isDirected():
			classConnectedComponents = properties.StronglyConnectedComponents
		else:
			classConnectedComponents = properties.ConnectedComponents

		for parameter in [
			("Node Centrality",	"Degree",						True,	funcScores,	"Score",				centrality.DegreeCentrality, 			(G, )),
			("Node Centrality",	"K-Core Decomposition",			True,	funcScores,	"Score",				centrality.CoreDecomposition, 			(G, )),
			("Node Centrality",	"Local Clustering Coefficient",	True,	funcScores,	"Score",				centrality.LocalClusteringCoefficient,	(G, )),
			("Node Centrality",	"PageRank",					True,	funcScores,	"Score",				centrality.PageRank, 					(G, )),
			("Node Centrality",	"K-Path Centrality",			True,	funcScores,	"Score",				centrality.KPathCentrality,				(G, )),
			("Node Centrality",	"Katz Centrality",				True,	funcScores,	"Score",				centrality.KatzCentrality,				(G, )),
			("Node Centrality",	"Betweenness (ap.)",					True,	funcScores,	"Score",				centrality.ApproxBetweenness2,			(G, max(42, math.log(G.numberOfNodes()), True))),
			("Node Centrality",	"Closeness",					True,	funcScores,	"Score",				centrality.ApproxCloseness,			(G, max(42, math.log(G.numberOfNodes()), True))),
			("Partition",		"Communities",					False,	funcSizes,	"Nodes per Community",	community.PLM,			 				(G, )),
			("Partition",		"Connected Components",			False,	funcSizes,	"Nodes per Component",	classConnectedComponents,				(G, )),
			("Partition",		"K-Core Decomposition",			False,	funcSizes,	"Nodes per Shell",		centrality.CoreDecomposition, 			(G, ))
		]: result.__addMeasure(parameter, exclude)

		if cls.__verbose:
			timerAll = stopwatch.Timer()
		result.__loadProperties()
		result.__loadMeasures()
		if cls.__verbose:
			if cls.__verboseLevel < 1:
				print("")
			print("\ntotal time (measures + stats + correlations): {:.2F} s".format(timerAll.elapsed))
			print("total speed: {:.1F} edges/s".format(G.numberOfEdges() / timerAll.elapsed))
		return result;


	@classmethod
	def setVerbose(cls, verbose=False, level=0, filename=""):
		""" TODO: """
		cls.__verbose = verbose
		cls.__verboseLevel = level
		cls.__verboseFilename = filename


	@classmethod
	def getVerbose(cls):
		""" TODO: """
		return (cls.__verbose, cls.__verboseLevel, cls.__verboseFilename)


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


	def output(self, type, directory, filename=None, style="light", color=(0, 0, 1), parallel=False):
		""" TODO """
		options_type = ["HTML", "LaTeX", None]
		for o in options_type:
			if o is None:
				raise ValueError("unknown type: options are " + str(options_type[0:len(options_type)-1]))
			if o == type:
				break;

		if type == "LaTeX":
			raise RuntimeError("not implemented, yet")

		result = self.__format(
			type = type,
			directory = directory,
			style = style,
			color = color,
			pageIndex = 0,
			parallel = parallel
		)

		if type == "HTML":
			js = readfile("html/profiling.js", False).replace("\\\\", "\\")
			css = readfile("html/profiling.css", False).replace("\\\\", "\\")
			result = """
				<!DOCTYPE HTML>
				<html>
					<head>
						<meta charset="utf-8">
						<style>
							""" + css + """
						</style>
						<script type="text/x-mathjax-config">
							MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$']]}});
						</script>
						<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
						<script>
							""" + js + """
						</script>
					</head>
					<body>
						""" + result + """
					</body>
				</html>
			"""
			if filename is None:
				filename  = "profile-{0}.html".format(self.__G.getName())
		else:
			raise Error("unknown output type")

		with open(directory + "/" + filename, 'w') as file:
			file.write(result)

	def show(self, style="light", color=(0, 0, 1), parallel=False):
		""" TODO: """
		try:
			__IPYTHON__
		except:
			raise RuntimeError("this function cannot be used outside ipython notebook")

		result = self.__format(
			type = "HTML",
			directory = "",
			style = style,
			color = color,
			pageIndex = self.__pageCount,
			parallel = parallel
		)

		display_html(HTML(result))
		self.__pageCount = self.__pageCount + 1


	def __format(self, type, directory, style, color, pageIndex, parallel):
		""" TODO """
		theme = plot.Theme()
		theme.set(style, color)

		pool = multiprocessing.ThreadPool(self.__parallel, parallel)
		for name in self.__measures:
			category = self.__measures[name]["category"]
			pool.put(
				plot.Measure(name, (
					0,
					self.__measures[name]["stat"],
					self.__measures[name]["label"],
					theme
				))
			)
			pool.put(
				plot.Measure(name, (
					1,
					self.__measures[name]["stat"],
					self.__measures[name]["label"],
					theme
				))
			)
			if category == "Partition":
				pool.put(
					plot.Measure(name, (
						2,
						self.__measures[name]["stat"],
						self.__measures[name]["label"],
						theme
					))
				)
		while pool.numberOfTasks() > 0:
			(plotType, name, data) = pool.get()
			try:
				category = self.__measures[name]["category"]

				if plotType == "Plot.Measure":
					(index, image) = data
					self.__measures[name]["image"][index] = image
			except Exception as e:
				self.__verbosePrint("Error (Post Processing): " + plotType + " - " + name, level=-1)
				self.___verbosePrint(str(e), level=-1)
		pool.join()

		if type == "HTML":
			templateProfile = readfile("html/profile.html", False)
			templateMeasure = readfile("html/measure.html", False)

		results = {}
		for category in self.__correlations:
			results[category] = {}
			results[category]["Correlations"] = {}
			results[category]["Correlations"]["HeatMaps"] = ""
			results[category]["Correlations"]["ScatterPlots"] = ""
			results[category]["Measures"] = ""
			results[category]["Overview"] = ""

			def funcHeatMap(category, correlationName):
				result = ""

				if type == "HTML":
					result += "<div class=\"SubCategory HeatTable\" data-title=\"" + correlationName + "\">"
					keyBList = []
					for keyA in self.__measures:
						if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
							keyBList.append(keyA)
							for keyB in keyBList:
								try:
									value = self.__correlations[category][keyA][keyB]
								except:
									value = self.__correlations[category][keyB][keyA]
								result += "<div class=\"HeatCell\" title=\"" + keyB + " - " + keyA + "\" data-heat=\"{:+.3F}\"></div>".format(value["stat"][correlationName])
							result += "<div class=\"HeatCellName\">" + keyB + "</div><br>"
					result += "</div>"

				return result
			results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Pearson's Correlation Coefficient")
			results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Spearman's Rank Correlation Coefficient")
			results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Fechner's Correlation Coefficient")

			def funcScatterPlot(category):
				result = ""

				if type == "HTML":
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

		# TODO: documentation
		for key in self.__measures:
			measure = self.__measures[key]
			name = measure["name"]
			category = measure["category"]
			image = measure["image"]
			stat = measure["stat"]
			assortativity = measure["assortativity"]
			centralization = measure["centralization"]


			description = "N/A"
			try:
				description = readfile("description/" + key + ".txt")
			except:
				pass

			extentions = ""
			try:
				if type == "HTML":
					extentions = "<div class=\"PartitionPie\"><img src=\"data:image/svg+xml;utf8," + image[2] + "\" /></div>"
			except:
				pass

			results[category]["Measures"] += self.__formatMeasureTemplate(
				templateMeasure,
				pageIndex,
				key,
				name,
				image,
				stat,
				assortativity,
				centralization,
				extentions,
				description
			)
			if type == "HTML":
				results[category]["Overview"] += "<div class=\"Thumbnail_Overview\" data-title=\"" + name + "\"><a href=\"#NetworKit_Page_" + str(pageIndex) + "_" + key + "\"><img src=\"data:image/svg+xml;utf8," + image[1] + "\" /></a></div>"

		result = self.__formatProfileTemplate(
			templateProfile,
			pageIndex,
			results
		)
		return result


	def __formatMeasureTemplate(self, template, pageIndex, key, name, image, stat, assortativity, centralization, extentions, description):
		""" TODO: """
		result = template.format(**locals())
		return result


	def __formatProfileTemplate(self, template, pageIndex, results):
		""" TODO: """
		properties = self.__properties
		result = template.format(**locals())
		return result


	def __addMeasure(self, args, exclude):
		""" TODO: """
		(measureCategory, measureName, correlate, getter, label, measureClass, parameters) = args
		measureKey = measureClass.__name__
		if measureKey not in exclude:
			measure = {}
			measure["name"] = measureName
			measure["category"] = measureCategory
			measure["correlate"] = correlate
			measure["getter"] = getter
			measure["label"] = label
			measure["class"] = measureClass
			measure["parameters"] = parameters
			measure["data"] = {}
			measure["image"] = {}
			self.__measures[measureKey] = measure
		try:
			self.__correlations[measureCategory]
		except:
			self.__correlations[measureCategory] = {}


	def __loadProperties(self):
		""" TODO: """
		self.__properties["Name"] = self.__G.getName()
		self.__properties["Nodes"] = self.__G.numberOfNodes()
		self.__properties["Edges"] = self.__G.numberOfEdges()
		self.__properties["Density"] = self.__G.density()
		self.__properties["Directed"] = self.__G.isDirected()
		self.__properties["Weighted"] = self.__G.isWeighted()
		self.__properties["Self Loops"] = self.__G.numberOfSelfLoops()

		timerInstance = stopwatch.Timer()
		self.__verbosePrint("Diameter: ", end="")
		try:
			diameter = properties.Diameter.estimatedDiameterRange(self.__G, error=0.1)
		except:
			diameter = "N/A"
		elapsedMain = timerInstance.elapsed
		self.__verbosePrint("{:.2F} s".format(elapsedMain))
		self.__verbosePrint("")
		self.__properties["Diameter Range"] = diameter


	def __loadMeasures(self):
		""" TODO: """
		pool = multiprocessing.ThreadPool(self.__parallel, False)

		for name in self.__measures:
			measure = self.__measures[name]
			self.__verbosePrint(name + ": ", end="")
			try:
				instance = measure["class"](*measure["parameters"])
			except Exception as e:
				del self.__measures[name]
				self.__verbosePrint("(removed)\n>> " + str(e))
				continue

			timerInstance = stopwatch.Timer()
			instance.run()
			measure["data"]["sample"] = measure["getter"](instance)
			elapsedMain = timerInstance.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedMain))

			self.__verbosePrint("    Sort: ", end="")
			timerPostSort = stopwatch.Timer()
			measure["data"]["sorted"] = stat.sorted(measure["data"]["sample"])
			elapsedPostSort = timerPostSort.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedPostSort))

			self.__verbosePrint("    Rank: ", end="")
			timerPostRank = stopwatch.Timer()
			measure["data"]["ranked"] = stat.ranked(measure["data"]["sample"])
			elapsedPostRank = timerPostRank.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedPostRank))

			self.__verbosePrint("    Assortativity: ", end="")
			timerPostAssortativity = stopwatch.Timer()
			if self.__measures[name]["category"] == "Node Centrality":
				assortativity = properties.Assortativity(self.__G, measure["data"]["sample"])
				assortativity.run()
				measure["assortativity"] = assortativity.getCoefficient()
			else:
				measure["assortativity"] = float("nan")
			elapsedPostAssortativity = timerPostAssortativity.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedPostAssortativity))

			# self.__verbosePrint("    Centralization: ", end="")
			timerPostCentralization = stopwatch.Timer()
			if self.__measures[name]["category"] == "Node Centrality":
				try:
					measure["centralization"] = instance.centralization()
				except:
					self.__verbosePrint("Centrality.centralization not properly defined for {0}".format(name), level=0)
					measure["centralization"] = float("nan")
			else:
				measure["centralization"] = float("nan")
			elapsedPostCentralization = timerPostCentralization.elapsed

			measure["time"] = (
				elapsedMain,
				elapsedPostSort,
				elapsedPostRank,
				elapsedPostAssortativity,
			)

		self.__verbosePrint("")

		for name in self.__measures:
			if len(self.__measures[name]["data"]["sample"]) <= 1:
				del self.__measures[name]
			else:
				category = self.__measures[name]["category"]
				pool.put(
					stat.Stat(name, (
						self.__measures[name]["data"]["sample"],
						self.__measures[name]["data"]["sorted"],
						self.__measures[name]["data"]["ranked"],
						category == "Partition"
					))
				)

		while pool.numberOfTasks() > 0:
			(type, name, data) = pool.get()

			try:
				category = self.__measures[name]["category"]

				if type == "Stat":
					self.__measures[name]["stat"] = data
					self.__verbosePrint("Stat: " + name, level=1)
					if self.__measures[name]["correlate"]:
						for key in self.__correlations[category]:
							self.__correlations[category][key][name] = {}
							self.__correlations[category][key][name]["stat"] = {}
							pool.put(
								stat.Correlation(key, (
									name,
									self.__measures[key]["data"]["sample"],
									self.__measures[key]["data"]["ranked"],
									self.__measures[key]["stat"],
									self.__measures[name]["data"]["sample"],
									self.__measures[name]["data"]["ranked"],
									self.__measures[name]["stat"]
								))
							)
							pool.put(
								plot.Scatter(key, (
									name,
									self.__measures[key]["label"],
									self.__measures[name]["label"],
									self.__measures[key]["data"]["sample"],
									self.__measures[name]["data"]["sample"]
								))
							)
						self.__correlations[category][name] = {}
						self.__correlations[category][name][name] = {}
						self.__correlations[category][name][name]["stat"] = {
							"Spearman's Rank Correlation Coefficient": 1,
							"Pearson's Correlation Coefficient": 1,
							"Fechner's Correlation Coefficient": 1
						}
						self.__correlations[category][name][name]["image"] = ""

				elif type == "Correlation":
					(nameB, correlation) = data
					self.__verbosePrint("Correlation: " + name + " <-> " + nameB, level=1)
					self.__correlations[category][name][nameB]["stat"] = correlation

				elif type == "Plot.Scatter":
					(nameB, image) = data
					self.__verbosePrint("Plot.Scatter: " + name, level=1)
					self.__correlations[category][name][nameB]["image"] = image
			except Exception as e:
				self.__verbosePrint("Error (Post Processing): " + type + " - " + name, level=-1)
				self.__verbosePrint(str(e), level=-1)

		pool.join()


	def __verbosePrint(self, text, end="\n", level=0):
		if self.__verboseLevel >= level:
			text = text + end
		else:
			text = "."

		if self.__verbose or level < 0:
			print(text, end="", flush=True)

		if self.__verboseFilename != "":
			with open(self.__verboseFilename, 'a+') as file:
				file.write(text)
