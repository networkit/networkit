from networkit import *
# import tabulate
# import pandas
# import seaborn
# import matplotlib.pyplot as plt
# from matplotlib._pylab_helpers import Gcf
# from IPython.core.pylabtools import print_figure
from IPython.core.display import *
from urllib.parse import quote
import io
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
		
		
def readfile(postfix):
	with open(__file__[:__file__.rfind(".py")] + "." + postfix, "r") as file:
		return " ".join(file.read().split())
		

def header(tag, type, data):
	result = """
		if (!document.getElementById('NetworKit_""" + tag + """')) {
			var element = document.createElement('""" + tag + """');
			element.type = 'text/""" + type + """';
			element.innerHTML = '""" + data + """';
			element.setAttribute('id', 'NetworKit_""" + tag + """');
			document.head.appendChild(element);
		} else {
			document.getElementById('NetworKit_""" + tag + """').innerHTML = '""" + data + """';
		}
	"""
	return result
	
	
display_html(
	HTML("""
		<script type="text/javascript">
		<!--
			""" + header("script", "javascript", readfile("js"))  + """
			""" + header("style",  "css",        readfile("css")) + """
		-->
		</script>
	""")
)


count = 0


def createImageURI():
	fig = plt.gcf()
	imgdata = io.StringIO()
	fig.savefig(imgdata, format='svg')
	plaintext = imgdata.getvalue()
	plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
	encoded = quote(plaintext, safe='');
	encoded = "data:image/svg+xml;utf8," + encoded; 
	return encoded;

	
def histogram(data, xlabel, ylabel):
	plt.clf();
	
	mu, sigma = 100, 15
	x = mu + sigma*np.random.randn(10000)

	# the histogram of the data
	n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
	
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	# plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
	plt.axis([40, 160, 0, 0.03])
	plt.grid(True)
	
	return createImageURI();

	
def plot2():
	plt.clf();
	
	plt.plot(range(30, 20))
	
	return createImageURI();

	
def profile(G):
	global count
	plt.ioff()
	
	degree = properties.degreeSequence(G)
	
	pageIndex = count
	centrality_1_images = histogram(degree, "x-Axis", "y-Axis");
	centrality_2_images = plot2();
	
	result = readfile("html")
	result = result.format(**locals());
	display_html(HTML(result))
	
	count = count + 1;
	
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
	
	
# def asSlideShow(*images):
	# """
	# """
	# global imageIndex
	# result = "<div>"
	# offset = imageIndex
	# for image in images:
		# result += "<div id='profiling_slideshow_img_" + str(imageIndex) + "'>" + image + "</div>"
		# imageIndex = imageIndex + 1
	# result += "<a href='javascript:SlideShow(" + str(offset) + "," + str(len(images)) + ");'><div>Button</div></a>"
	# result += "</div>"
	# return result 


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