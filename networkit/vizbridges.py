# local imports
from .support import MissingDependencyError
from .structures import Partition
from .viz import MaxentStress

# external imports
try:
	import ipycytoscape
except ImportError:
	hasCyto = False
else:
	hasCyto = True

try:
	import plotly.graph_objs as go
	import plotly.express as px
except ImportError:
	hasPlotly = False
else:
	hasPlotly = True
try:
	import seaborn
except ImportError:
	hasSeaborn = False
else:
	hasSeaborn = True

class Dimension:
	"""
	Supported dimensions for visualization.

	Possible values:

	- networkit.vizbridges.Dimension.Two (visualization in 2D with Cytoscape)
	- networkit.vizbridges.Dimension.TwoForcePlotly (visualization in 2D with Plotly)
	- networkit.vizbridges.Dimension.Three (visualization in 3D with Plotly)
	"""
	Two = 0
	TwoForcePlotly = 1
	Three = 2

def _getColorPalette(nodePalette = None, nodePartition = None):
	# Set color palettes
	if nodePalette is not None:
		palette = nodePalette
	else:
		if not hasSeaborn:
			raise MissingDependencyError("seaborn")
		# Partitions and scores have different default color palettes
		if nodePartition is not None:
			palette = seaborn.color_palette("hls", nodePartition.numberOfSubsets())
		else:
			palette = seaborn.color_palette("rocket_r", as_cmap=True).colors
	return palette

def _calculateNodeColoring(G, palette, nodeScores = None, nodePartition = None):

	# Color calculation: score = continuous distribution, partition = discrete distribution
	hcColors = []

	# Partition
	if nodePartition is not None:
		if len(palette) < nodePartition.numberOfSubsets():
			raise Exception("Number of partitions higher than number of colors in provided palette. Provide node_palette with enough colors.")

		partitions = nodePartition.getVector()

		if len(palette) < nodePartition.numberOfSubsets():
			raise Exception("Number of partitions to high for default coloring. Provide node_palette with enough colors.")

		for i in range(0, len(partitions)):
			hcColors.append((palette[partitions[i]][0] * 255, palette[partitions[i]][1] * 255, palette[partitions[i]][2] * 255))

	# Score
	elif nodeScores is not None:

		minhc = min(nodeScores)
		maxhc = max(nodeScores)

		# Calculate coloring of nodes
		def getRgb(minimum, maximum, value):
			minimum, maximum, value = float(minimum), float(maximum), float(value)
			ratio = int((len(palette) - 1) * (value-minimum) / (maximum - minimum) * (value-minimum) / (maximum - minimum))
			r = int(palette[ratio][0] * 255)
			g = int(palette[ratio][1] * 255)
			b = int(palette[ratio][2] * 255)
			return r, g, b

		if abs(maxhc - minhc) > 0:
			for score in nodeScores:
				hcColors.append(getRgb(minhc, maxhc, score))
		else:
			color = palette[int((len(palette) -1) / 2)];
			for i in range(0, len(nodeScores)):
				hcColors.append((color[0] * 255, color[1] * 255, color[2] * 255))

	# No node values
	else:
		color = palette[0];
		for i in G.iterNodes():
			hcColors.append((color[0] * 255, color[1] * 255, color[2] * 255))

	return hcColors


def widgetFromGraph(G, dimension = Dimension.Two, nodeScores = None, nodePartition = None, nodePalette = None, showIds = True, customSize = None):
	""" 
	widgetFromGraph(G, dimension=Dimension.Two, nodeScores=None, nodePartition=None, nodePalette=None, showIds=True, customSize=None)

	Creates a widget with a visualization of a given graph. The widget uses one of the supported
	plugins - either Cytoscape (2D) or Plotly (3D). The returned widget already contains
	all nodes and edges from the graph. The graph is colored using an array of norm. rgb-values
	based on seaborn perceptually uniform color map (rocket) or a user given custom color array. 
	See matplotlib color maps for the correct formatting: 
	https://matplotlib.org/api/_as_gen/matplotlib.colors.Colormap.html#matplotlib.colors.Colormap


	Parameters
	----------
	G : networkit.graph.Graph
		The input graph.
	dimension : networkit.vizbridges.Dimension, optional
		Select whether to plot in 2D or 3D. This also influences which plugin is choosen
		for visualization. For 2D Cytoscape with auto-layouting is used, for 3D Plotly with
		a Maxent-Stress layout is used. Option :code:`Dimension.TwoForcePlotly` forces a plot 
		in 2D with using Plotly instead of Cytoscape. Default: networkit.vizbridges.Dimension.Two
	nodeScores : list(float), optional
		List of scores for each nodes, for example scores from a centrality measure. This is 
		used for color-calculation of the nodes (continuous distribution). Provide either 
		node_scores or node_partition - not both. Default: None
	nodePartition : networkit.structures.Partition, optional
		Partition object. This is used for color-calculation of the nodes (discrete distribution). 
		Provide either node_scores or node_partition - not both. Default: None	
	nodePalette : list(tuple(float, float, float)), optional
		List consisting of normalized rgb-values. If none is given, seaborn.color_palette.colors is used. Default: None
	showIds : boolean, optional
		Set whether node ids should be visible in plot-widget. Default: True
	customSize : int, optional
		If not set, plugins will use a default size. Otherwise the widget will have a certain width and height. Default: None       
	"""
	# Sanity checks
	if dimension == Dimension.Two and not hasCyto:
		raise MissingDependencyError("ipycytoscape")

	if (dimension == Dimension.Three or dimension == Dimension.TwoForcePlotly) and not hasPlotly:
		raise MissingDependencyError("Plotly")

	if nodeScores is not None:
		if nodePartition is not None:
			raise Exception("Provide either nodeScores or nodePartition - not both.")
		if len(nodeScores) != G.upperNodeIdBound():
			raise Exception("nodeScores should include scores for every node.")	

	# Color palette is needed for node coloring with Plotly and Cytoscape
	palette = _getColorPalette(nodePalette=nodePalette, nodePartition=nodePartition)
	
	# Set styling
	if showIds:
		s = [{
			'selector': 'node',
			'css': {
				'background-color': 'data(color)',
				'content': 'data(id)'}
			}]
	else:
		s = [{
			'selector': 'node',
			'css': {
				'background-color': 'data(color)'}
			}]

	if dimension == Dimension.Two:
		# Create widget
		graphWidget = ipycytoscape.CytoscapeWidget()

		# Add data
		nodes = []
		edges = []
		if G.isDirected():
			edgeClass = "directed "
		else:
			edgeClass = "undirected "

		# Color list (norm. RGB-values) is needed for node coloring with Cytoscape 
		hcColors = _calculateNodeColoring(G, palette, nodeScores, nodePartition)

		for i in G.iterNodes():
			n = ipycytoscape.Node(data={"id": i, "color": hcColors[i] })						
			nodes.append(n)

		for u,v in G.iterEdges():
			e = ipycytoscape.Edge(data={"source": u, "target": v, "classes": edgeClass })
			edges.append(e)

		# It is much faster to add edges and nodes in bulk.
		graphWidget.graph.add_nodes(nodes)
		graphWidget.graph.add_edges(edges, G.isDirected())

		# Set layout
		graphWidget.set_style(s)
		graphWidget.set_layout(name='cose')

	else:
		# Create widget
		graphWidget = go.FigureWidget()

		# Set layout
		maxLayout = MaxentStress(G, 3, 3, fastComputation=1, graphDistance=0)
		maxLayout.run()
		coordinates = maxLayout.getCoordinates()

		# Set node coloring
		if nodePartition is not None:
			scores = nodePartition.getVector()
		elif nodeScores is not None:
			scores = nodeScores
		else:
			scores = [0.0] * G.upperNodeIdBound()
		labels = ["Node: " + str(id) + "<br>Score: " + str(score) for id, score in enumerate(scores)]

		# Initiate widget data
		nodes = [[],[],[]]
		nodes[0] = [coordinates[k][0] for k in G.iterNodes()]
		nodes[1] = [coordinates[k][1] for k in G.iterNodes()]
		nodes[2] = [coordinates[k][2] for k in G.iterNodes()]

		index = 0
		if dimension == Dimension.TwoForcePlotly:
			edges = [[None] * G.numberOfEdges() * 2,[None] * G.numberOfEdges() * 2, [None] * G.numberOfEdges() * 2]
			for e in G.iterEdges():
				edges[0][index] = coordinates[e[0]][0]
				edges[0][index+1] = coordinates[e[1]][0]
				edges[1][index] = coordinates[e[0]][1]
				edges[1][index+1] = coordinates[e[1]][1]
				index = index + 2
		
			nodeScatter = go.Scatter(x=nodes[0],
				y=nodes[1],
				mode='markers',
				name='nodes',
				marker=dict(symbol='circle',
					size=9,
					colorscale=px.colors.convert_colorscale_to_rgb(px.colors.make_colorscale(palette)),
					color=scores,
					line=dict(color='rgb(50,50,50)', width=0.5)),
				hoverinfo='text',
				text = labels)

			edgeScatter = go.Scatter(x=edges[0],
				y=edges[1],
				mode='lines',
				opacity=0.7,
				line= dict(color='rgb(180,180,180)', width=2),
				hoverinfo='none',
				showlegend=None,
				name='edges')

		else:
			edges = [[None] * G.numberOfEdges() * 3,[None] * G.numberOfEdges() * 3, [None] * G.numberOfEdges() * 3]
			for e in G.iterEdges():
				edges[0][index] = coordinates[e[0]][0]
				edges[0][index+1] = coordinates[e[1]][0]
				edges[1][index] = coordinates[e[0]][1]
				edges[1][index+1] = coordinates[e[1]][1]
				edges[2][index] = coordinates[e[0]][2]
				edges[2][index+1] = coordinates[e[1]][2]
				index = index + 3

			nodeScatter = go.Scatter3d(x=nodes[0],
				y=nodes[1],
				z=nodes[2],
				mode='markers',
				name='nodes',
				marker=dict(symbol='circle',
					size=9,
					colorscale=px.colors.convert_colorscale_to_rgb(px.colors.make_colorscale(palette)),
					color=scores,
					line=dict(color='rgb(50,50,50)', width=0.5)),
				hoverinfo='text',
				text = labels)

			edgeScatter = go.Scatter3d(x=edges[0],
				y=edges[1],
				z=edges[2],
				mode='lines',
				opacity=0.7,
				line= dict(color='rgb(180,180,180)', width=2),
				hoverinfo='none',
				showlegend=None,
				name='edges')

		graphWidget.add_traces(nodeScatter)
		graphWidget.add_traces(edgeScatter)       

		# Set layout
		minCoordinate, maxCoordinate = 0.0, 0.0
		for pos in coordinates:
			minCoordinate = min(minCoordinate, min(pos))
			maxCoordinate = max(maxCoordinate, max(pos))

		if customSize:
			width, height = customSize, customSize
		else:
			width, height = 1000, 1000

		axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
			zeroline=False,
			showgrid=True,
			showticklabels=True,
			title='',
			autorange=False,
			range=[1.1 * minCoordinate, 1.1 * maxCoordinate])

		graphWidget.layout = go.Layout(font= dict(size=12),
			showlegend=False,
			autosize=True,
			scene_aspectmode='cube',
			width=width,
			height=height,
			scene=dict(
				xaxis=dict(axis),
				yaxis=dict(axis),
				zaxis=dict(axis),
				bgcolor='white',
				camera=dict(
					eye={'x' : 1.3, 'y' : 1.3, 'z' : 0.5})
			),
			margin=go.layout.Margin(
				l=10,
				r=10,
				b=10,
				t=10,
			),
			hovermode='closest')
	return graphWidget
