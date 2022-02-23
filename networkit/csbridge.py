# local imports
from .support import MissingDependencyError
from .structures import Partition
from .vizbridges import widgetFromGraph

# external imports
try:
	import ipycytoscape
except ImportError:
	have_cyto = False
else:
	have_cyto = True

try:
	import seaborn
except ImportError:
	have_seaborn = False
else:
	have_seaborn = True

def widget_from_graph(G, node_scores = None, node_partition = None, node_palette = None, show_ids = True):
	""" 
	widget_from_graph(G, node_scores = None, node_partition = None, node_palette = None, show_ids = True)

	Note
	----
	DEPRECATED. Use networkit.vizbridges.widgetFromGraph() instead.

	Creates a ipycytoscape-widget from a given graph. The returned widget already contains
	all nodes and edges from the graph. The graph is colored using an array of norm. rgb-values
	based on seaborn perceptually uniform color map (rocket) or a user given custom color array. 
	See matplotlib color maps for the correct formatting: 
	https://matplotlib.org/api/_as_gen/matplotlib.colors.Colormap.html#matplotlib.colors.Colormap

	Parameters
	----------
	G : networkit.graph.Graph
		Graph (nodes, edges), which should be visualized
	node_scores : list of numbers, optional
		List of scores for each nodes, for example scores from a centrality measure. This is 
		used for color-calculation of the nodes (continuous distribution). Provide either 
		node_scores or node_partition - not both. Default: None
	node_partition : networkit.structures.Partition, optional
		Partition object. This is used for color-calculation of the nodes (discrete distribution). 
		Provide either node_scores or node_partition - not both. Default: None 	
	node_palette : list of tuples, optional
		Array consisting of normalized rgb-values. If none is given, seaborn.color_palette.colors is used.  Default: None
	show_ids : 	boolean, optional
		Set whether node ids should be visible in plot-widget. Is set to True by default. Default: None
	"""
	print("WARNING: Module csbridge is deprecated and will be removed in future updates. Use networkit.vizbridges.widgetFromGraph() instead.")
	return widgetFromGraph(G, nodeScores = node_scores, nodePartition = node_partition, nodePalette = node_palette, showIds = show_ids)
