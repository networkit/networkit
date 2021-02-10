# local imports
from .support import MissingDependencyError
from .structures import Partition

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

def widget_from_graph(G, node_scores = None, node_partition = None, node_colors = None, show_ids = True):
	""" 
	Creates a ipycytoscape-widget from a given graph. The returned widget already contains
	all nodes and edges from the graph. The graph is colored using an array of norm. rgb-values
	based on seaborn perceptually uniform color map (rocket) or a user given custom color array. 
	See matplotlib color maps for the correct formatting: 
	https://matplotlib.org/api/_as_gen/matplotlib.colors.Colormap.html#matplotlib.colors.Colormap 

	Parameters:
	-----------
	G : networkit.graph.Graph
		Graph (nodes, edges), which should be visualized
	node_scores : list of numbers
		List of scores for each nodes, for example scores from a centrality measure. This is 
		used for color-calculation of the nodes (continuous distribution). Provide either 
		node_scores or node_partition - not both.
	node_partition : networkit.structures.Partition
		Partition object. This is used for color-calculation of the nodes (discrete distribution). 
		Provide either node_scores or node_partition - not both. 	
	node_colors : list of tuples
		Array consisting of normalized rgb-values. If none is given, seaborn.color_palette.colors is used
	show_ids : 	boolean	
		Set whether node ids should be visible in plot-widget. Is set to True by default. 
	"""
	# Sanity checks
	if not have_cyto:
		raise MissingDependencyError("ipycytoscape")

	if node_scores is not None and node_partition is not None:
		raise InputError("Provide either node_scores or node_partition - not both")

	# Color calculation: score = continuous distribution, partition = discrete distribution
	hc_colors = []

	# Partition
	if node_partition is not None:
		if node_colors is not None:
			palette = node_colors
		else:
			if not have_seaborn:
				raise MissingDependencyError("seaborn")
			palette = seaborn.color_palette("hls", node_partition.numberOfSubsets())

		partitions = node_partition.getVector()

		for i in G.iterNodes():
			hc_colors.append((palette[partitions[i]][0] * 255, palette[partitions[i]][1] * 255, palette[partitions[i]][2] * 255))

	# Score
	elif node_scores is not None and len(node_scores) == G.numberOfNodes():
		if node_colors is not None:
			palette = node_colors
		else:
			if not have_seaborn:
				raise MissingDependencyError("seaborn")
			palette = seaborn.color_palette("rocket_r", as_cmap=True).colors

		minhc = min(node_scores)
		maxhc = max(node_scores)

		# calculate coloring of nodes
		def get_rgb(minimum, maximum, value):
			minimum, maximum, value = float(minimum), float(maximum), float(value)
			ratio = int((len(palette) - 1) * (value-minimum) / (maximum - minimum) * (value-minimum) / (maximum - minimum))
			r = int(palette[ratio][0] * 255)
			g = int(palette[ratio][1] * 255)
			b = int(palette[ratio][2] * 255)
			return r, g, b

		if abs(maxhc - minhc) > 0:
			for score in node_scores:
				hc_colors.append(get_rgb(minhc, maxhc, score))
		else: 
			for i in G.iterNodes():
				hc_colors.append((247, 137, 97))

	# No node values
	else:
		for i in G.iterNodes():
			hc_colors.append((247, 137, 97))
	
	# Set styling
	if show_ids:
		s = [{
			'selector': 'node',
				'css': {
					'background-color': 'data(color)',
					'content': 'data(id)'
				}}]
	else:
		s = [{
			'selector': 'node',
				'css': {
					'background-color': 'data(color)'
				}}]

	# Create widget
	cytoWidget = ipycytoscape.CytoscapeWidget()

	nodes = []
	edges = []

	if G.isDirected():
		edge_class = "directed "
	else:
		edge_class = "undirected "

	for i in G.iterNodes():
		n = ipycytoscape.Node(data={"id": i, "color": hc_colors[i] })
		nodes.append(n)

	for u,v in G.iterEdges():
		e = ipycytoscape.Edge(data={"source": u, "target": v, "classes": edge_class })
		edges.append(e)

	# It is much faster to add edges and nodes in bulk.
	cytoWidget.graph.add_nodes(nodes)
	cytoWidget.graph.add_edges(edges, G.isDirected())

	cytoWidget.set_style(s)
	cytoWidget.set_layout(name='cose')
	return cytoWidget
