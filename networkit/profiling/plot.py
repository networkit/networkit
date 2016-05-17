#
# file: plot.py
# author: Mark Erb
#

from . import job

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages

import io
from urllib.parse import quote
import math


class Theme:
	""" layout theme for plots """
	
	def __init__(self):
		""" constructor """
		self.set()

	@classmethod
	def RGBA2RGB(cls, color, alpha, background):
		""" converts a color with an given alpha to RGB for an fixed RGB background color """
		result = (
			color[0] * alpha + background[0] * (1-alpha),
			color[1] * alpha + background[1] * (1-alpha),
			color[2] * alpha + background[2] * (1-alpha),
			1
		)
		return result


	def set(self, style="light", color=(0, 0, 1)):
		""" sets style and color of the theme 
		
		Args:
			style: ("light")
			color: RGB tuple
		"""
		optionsStyle = ["light", "system"]
		if style not in optionsStyle:
			raise ValueError("possible style options: " + str(optionsStyle))
		if len(color) != 3:
			raise ValueError("(r,g,b) tuple required")

		if style == "system":
			self.__rcParams = mpl.rcParams
			raise ValueError("not implemented, yet")

		if style == "light":
			self.__defaultColor = (0, 0, 0)
			self.__defaultWidth = 1
			self.__backgroundColor = (1, 1, 1)
			self.__plotColor = Theme.RGBA2RGB(color, 0.6, self.__backgroundColor)
			self.__plotWidth = 3
			self.__faceColor = (color[0], color[1], color[2], 0.2)
			self.__faceColorGray = "lightgray"
			self.__edgeColor = (color[0], color[1], color[2], 0.6)
			self.__edgeColorGray = (0, 0, 0)
			self.__edgeWidth = 1
			self.__gridColor = "lightgray"
			self.__fontColor = (0, 0, 0)
			self.__fontSize = 10

		self.__color = color
		self.__style = style


	def get(self):
		""" return style and color """
		return (self.__style, self.__color)

		
	def getRcParams():
		""" return  matlibplot system parameters used in the theme """
		return self.__rcParams

		
	def getDefaultColor(self):
		""" returns the default color value of the theme """
		return self.__defaultColor
		
		
	def getDefaultWidth(self):
		""" returns the default width value of the theme """
		return self.__defaultWidth
		
		
	def getPlotColor(self):
		""" returns the plot color value of the theme """
		return self.__plotColor
		
		
	def getPlotWidth(self):
		""" returns the plot width value of the theme """
		return self.__plotWidth
		
		
	def getFaceColor(self):
		""" returns the face color value of the theme """
		return self.__faceColor
		
		
	def getFaceColorGray(self):
		""" returns the face color (gray) value of the theme """
		return self.__faceColorGray
		
		
	def getEdgeColor(self):
		""" returns the edge color value of the theme """
		return self.__edgeColor
		
		
	def getEdgeColorGray(self):
		""" returns the edge color (gray) value of the theme """
		return self.__edgeColorGray
		
		
	def getEdgeWidth(self):
		""" returns the edge width value of the theme """
		return self.__edgeWidth
		
		
	def getBackgroundColor(self):
		""" returns the background color value of the theme """
		return self.__backgroundColor
		
		
	def getGridColor(self):
		""" returns the grid color value of the theme """
		return self.__gridColor
		
		
	def getFontSize(self):
		""" returns the font size value of the theme """
		return self.__fontSize
		
		
	def getFontColor(self):
		""" returns the font color value of the theme """
		return self.__fontColor

		
class PlotJob(job.Job):
	def __init__(self, typename, plottype, options, name, params):
		""" constructor 
		
		Arg:
			typename: job type name
			plottype: "SVG" / "PDF"
			options: format dependent
			name: name of the measure
			params: a tuple of additional parameters
		"""
		job.Job.__init__(
			self,
			typename,
			name
		)
		self.__plottype = plottype
		self.__options = options
		self.__params = params
		
		
	def getParams(self):
		""" returns params """
		return self.__params		
	
	
	def save(self, id, fig, extention):
		""" generate plot output """
		result = ""

		fig.tight_layout()

		if self.__plottype == "SVG":
			imgdata = io.StringIO()
			fig.savefig(imgdata, format='svg')

			plaintext = imgdata.getvalue()
			plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
			result = quote(plaintext, safe='')

		elif self.__plottype == "PDF":
			filename = self.__options[1] + extention
			with PdfPages(self.__options[0] + "/" + filename + ".pdf") as pdf:
				pdf.savefig(fig)
			result = filename

		else:
			pass

		plt.close(fig)
		return (id, result)


		
class Measure(PlotJob):
	""" plot generation object for measures """
	
	def __init__(self, plottype, options, name, params):
		""" constructor: see PlotJob and .run() """
		PlotJob.__init__(
			self,
			"Plot.Measure",
			plottype,
			options,
			name,
			params
		)


	def run(self):
		""" computation """
		(index, stat, category, name, label, theme) = self.getParams()
		plt.ioff()

		plottype = "plot"

		def funcSpace(min, max):
			""" returns space within plot """
			result = 0.1
			if min < max:
				result = (max - min) * 0.04
			return result


		def funcTicks(min, max, numberOfTicks):
			""" returns ticks """
			result = []
			if numberOfTicks > 0:
				value = min
				step = (max - min) / numberOfTicks
				for i in range(numberOfTicks):
					result.append(min + step*i)
					value += step
				result.append(max)
			return result


		def funcPlotEnd(fig, ax, theme, width, height, x_showTicks=True, x_showTickLabels=True, y_showTicks=True, y_showTickLabels=True, drawAxis=True, showGrid=True):
			""" set some layout options """
			if not x_showTicks:
				ax.set_xticks([])
			if not x_showTickLabels:
				ax.set_xticklabels([])
			else:
				xfmt = ScalarFormatter(useMathText=True)
				xfmt.set_powerlimits((-2,3))
				ax.xaxis.set_major_formatter(xfmt)
			if not y_showTicks:
				ax.set_yticks([])
			if not y_showTickLabels:
				ax.set_yticklabels([])
			ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
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


		def funcPlotBox(ax):
			""" Box Plot """
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
			return ax


		def funcPlotPDF(ax):
			""" Histogram Plot """
			numberOfBins = stat["Binning"]["Number Histogram"]
			intervals = stat["Binning"]["Intervals Histogram"]
			absoluteFrequencies = stat["Binning"]["Absolute Frequencies Histogram"]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			y_min = 0
			y_max = stat["Binning"]["Mode"][1]
			x_space = funcSpace(x_min, x_max)
			y_space = funcSpace(y_min, y_max)
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
			return ax


		def funcPlotCDF(ax):
			""" Histogram Plot (commulativ) """
			numberOfBins = stat["Binning"]["Number CDF"]
			intervals = stat["Binning"]["Intervals CDF"]
			comulativeRelativeFrequencies = stat["Binning"]["Relative Frequencies CDF"]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			x_space = funcSpace(x_min, x_max)
			y_space = funcSpace(0, 1)
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
			return ax


		def funcPlotPie(ax):
			""" Pie Plot """
			numberOfTooSmallSubsets = stat["Binning"]["Pie"][1]
			relativeFrequencies = stat["Binning"]["Pie"][0]
			radius = 2
			accumulator = 0

			for i in range(len(relativeFrequencies)):
				if i == 0 and numberOfTooSmallSubsets == 0:
					continue

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
			ax1.set_ylabel('Box')
			ax1.yaxis.set_label_position("right")
			funcPlotBox(ax1)
			funcPlotEnd(
				fig = fig,
				ax = ax1,
				theme = theme,
				width = 4,
				height = 0.2,
				x_showTickLabels = False,
				y_showTicks = False,
				y_showTickLabels = False,
				showGrid = False
			)

			ax2 = plt.subplot2grid((40, 8), (3, 0), colspan=8, rowspan=20)
			ax2.set_ylabel("PDF (absolute)")
			ax2.yaxis.set_label_position("right")
			funcPlotPDF(ax2)
			funcPlotEnd(
				fig = fig,
				ax = ax2,
				theme = theme,
				width = 4,
				height = 3,
				x_showTickLabels = False
			)

			ax3 = plt.subplot2grid((40, 8), (23, 0), colspan=8, rowspan=17)
			ax3.set_xlabel(label)
			ax3.set_ylabel('CDF (relative)')
			ax3.yaxis.set_label_position("right")
			funcPlotCDF(ax3)
			funcPlotEnd(
				fig = fig,
				ax = ax3,
				theme = theme,
				width = 4,
				height = 3
			)
			fig.set_size_inches(6, 6)
			plottype = "stat"

		elif index == 1:
			fig = plt.figure()
			ax = fig.gca()

			ax.set_xlabel(label)
			# ax.xaxis.set_label_position("top")
			ax.set_ylabel("PDF (absolute)")
			ax.yaxis.set_label_position("right")
			funcPlotPDF(ax)
			funcPlotEnd(
				fig = fig,
				ax = ax,
				theme = theme,
				width = 4,
				height = 2.5
			)
			plottype = "thumb"

		elif index == 2:
			fig = plt.figure()
			ax = fig.gca()

			funcPlotPie(ax)
			funcPlotEnd(
				fig = fig,
				ax = ax,
				theme = theme,
				width = 6.4*1.5,
				height = 5.4*1.5,
				drawAxis = False
			)
			plottype = "pie"

		return self.save(
			index,
			fig,
			plottype + "-" + category + "-" + name
		)


class Scatter(PlotJob):
	""" plot generation object for scatter plots between measures """
	
	def __init__(self, plottype, options, name, params):
		""" constructor: see PlotJob and .run() """
		PlotJob.__init__(
			self,
			"Plot.Scatter",
			plottype,
			options,
			name,
			params
		)
		
		
	def run(self):
		""" computation """
		(name, nameA, nameB, labelA, labelB, stat_1, stat_2, correlation, theme) = self.getParams()
		plt.ioff()

		def funcHexBin(ax):
			gridsize = correlation["Binning"]["Grid Size"]
			frequencies = correlation["Binning"]["Absolute Frequencies"]
			max  = correlation["Binning"]["Max Frequency"]
			offsets = correlation["Binning"]["Offsets"]
			paths = correlation["Binning"]["Paths"]
			x_min = stat_1["Location"]["Min"]
			x_max = stat_1["Location"]["Max"]
			y_min = stat_2["Location"]["Min"]
			y_max = stat_2["Location"]["Max"]
			for i in range(len(frequencies)):
				color = Theme.RGBA2RGB(
					theme.getPlotColor(),
					math.log(frequencies[i]+1,10)/math.log(max+1,10),
					theme.getBackgroundColor()
				)
				path = paths.transformed(mpl.transforms.Affine2D().translate(
					offsets[i][0],
					offsets[i][1]
				)) 
				ax.add_patch(patches.PathPatch(
					path,
					facecolor = color,
					linestyle = "solid",
					linewidth = 0			
				))
			ax.set_xlim([x_min, x_max])
			ax.set_ylim([y_min, y_max])
			ax.set_xlabel(labelA)
			ax.set_ylabel(labelB)
			ax2 = ax.twinx()
			ax2.set_ylabel(nameB)
			ax2.set_yticks([])
			ax3 = ax.twiny()
			ax3.set_xlabel(nameA)
			ax3.set_xticks([])

		fig = plt.figure()
		ax = fig.gca()

		funcHexBin(ax)
		xfmt = ScalarFormatter(useMathText=True)
		xfmt.set_powerlimits((-1,1))
		ax.xaxis.set_major_formatter(xfmt)
		yfmt = ScalarFormatter(useMathText=True)
		yfmt.set_powerlimits((-1,1))
		ax.yaxis.set_major_formatter(yfmt)

		fig.set_size_inches(4, 3.75)

		return self.save(
			name,
			fig,
			"scatter." + nameA + " - " + nameB
		)