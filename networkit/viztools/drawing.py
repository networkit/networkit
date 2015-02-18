'''
Created on 09.04.2013
@author: cls

'''
# external imports
import networkx as nx
import matplotlib.pyplot as plt
import pylab


class GraphDrawer:
	
	def __init__(self, layout=nx.spring_layout, size=(8,8), nodeOpts={}, edgeOpts={}, labelOpts={}, labelled=False, emphasizeCliques=False):
		"""

		:param layout: The layout of the graph
		:param size: The size of the figure
		:param nodeOpts: additional node options
		:param edgeOpts: additional edge options
		:param labelOpts: additional label options
		:param labelled: node labels are drawn if true
		:param emphasizeCliques: emphasize Cliques of a graph
		"""
		self.layout = layout
		self.size = size
		self.labelled = labelled	# turn node labelling on or off
		self.emphCliques = emphasizeCliques


		# default options for drawing
		self.defaultEdgeOpts = {"edge_color": "darkgray"}
		self.defaultNodeOpts = {"node_color": "lightgray", "linewidths": 1.0, "node_size": (150 if self.labelled else 80)}
		self.defaultLabelOpts = {"font_size": 7}

		# if options are given as arguments, they will override the defaults
		self.edgeOpts = dict(self.defaultEdgeOpts, **edgeOpts)
		self.nodeOpts = dict(self.defaultNodeOpts, **nodeOpts)
		self.labelOpts = dict(self.defaultLabelOpts, **labelOpts)

		
		
		
	def draw(self, G, pos=None, nodeSizes=None, nodeColors=None):
		"""
		Draw graph unweighted

		:param G: The graph object
		:param pos: optional predefined node positions
		:param nodeSizes: optional array of node sizes
		"""
		
		plt.figure(figsize=self.size)
		plt.axis("off")
		if not pos:
			pos = self.layout(G)

		
		# draw nodes
		if nodeSizes:
			nodeOpts = dict(self.nodeOpts, **{"node_size" : nodeSizes})
		else:
			nodeOpts = self.nodeOpts

		if not nodeColors is None:
			nodeOpts = dict(nodeOpts, **{"node_color": nodeColors})

		nx.draw_networkx_nodes(G, pos, **nodeOpts)

		# draw nodelabels
		if self.labelled:
			nx.draw_networkx_labels(G, pos, **self.labelOpts)

		# draw edges
		if self.emphCliques:
			alpha = self.emphasizeCliques(G)
			i = 0
			for u,v in G.edges():
				nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], alpha=alpha[i])
				i += 1
		else:
			nx.draw_networkx_edges(G, pos, **self.edgeOpts)
		


	def emphasizeCliques(self, G, pos=None):
		""" Emphasize cliques of a graph.
		It creates a list with values for the edge transparency.
		The more the nodes have common neighbors, the darker will be the edges

		:param G: The graph object
		:param pos: optional predefined node positions
		:rtype: list
		"""

		plt.figure(figsize=self.size)
		plt.axis("off")
		if not pos:
			pos = self.layout(G)
			
		alpha = []

		# values for alpha are calculated depending on common neighbors between each two nodes
		for u,v in G.edges():
			nu = set(G.neighbors(u)+[u]) #neighbors from u
			nv = set(G.neighbors(v)+[v]) #neighbors from v
			intersection = nu.intersection(nv) #neighbors from both u and v
			union = nu.union(nv) #all neighbors
			average = (len(intersection)+len(union))/2
			alpha += [round(len(intersection)/average, 1)]

		# set minimum for edge transparency
		for i in range(len(G.edges())):
			if alpha[i] < 0.1:
				alpha[i] = 0.1

		return alpha


	def drawSelSCAN(self, G, seeds, cores, communities, pos=None):
		"""
		Draw graph unweighted and visualize seed and core nodes

		:param G: The graph object
		:param seeds: List of seed nodes
		:param cores: List of core nodes
		:param communities: Dictionary of communities. Nodes in the same community have the same color.
		:param pos: optional predefined node positions
		"""
		
		plt.figure(figsize=self.size)
		if not pos:
			pos = self.layout(G)
		plt.axis("off")

		colors = []
		for i in communities:
			colors += [communities[i]]
				
		# sort the nodes into different categories
		seedsAndCores = [u for u in G.nodes() if (u in seeds) and (u in cores)]
		onlySeeds = [u for u in G.nodes() if (u in seeds) and (u not in cores)]
		onlyCores = [u for u in G.nodes() if (u in cores) and (u not in seeds)]
		others = [u for u in G.nodes() if (u not in cores) and (u not in seeds)]
				
		edgeOpts = {"alpha":0.7}
		nodeOpts = {"node_color": [colors[u] for u in others], "alpha":0.8, "linewidths":3.5, "node_size": 350}
		corenodeOpts = {"node_color": [colors[u] for u in onlyCores], "alpha":0.8, "linewidths":2.5, "node_size": 420,"node_shape": 'p'}
		seednodeOpts = {"node_color":[colors[u] for u in onlySeeds],"alpha":0.8, "linewidths":2.5, "node_size": 350}
		seedcorenodeOpts = {"node_color":[colors[u] for u in  seedsAndCores], "alpha":0.8, "linewidths":3.5, "node_size": 420, "node_shape": 'p'}
		labelOpts = {"font_size": 10}

		# draw labels
		nx.draw_networkx_labels(G, pos, **labelOpts)
		
		# draw G.edges()
		coreEdges = [(u, v) for (u, v) in G.edges() if (u in cores) and (v in cores)]
		nonCoreEdges = [(u, v) for (u, v) in G.edges() if ((u, v) not in coreEdges)]	   
				
		nx.draw_networkx_edges(G, pos, edgelist=coreEdges, width=1.8)					#draw edges between cores
		nx.draw_networkx_edges(G, pos, edgelist=nonCoreEdges, **edgeOpts)
				

		
		# draw the nodes from the categories
		nx.draw_networkx_nodes(G, pos, nodelist=onlySeeds,**seednodeOpts)
		nx.draw_networkx_nodes(G, pos, nodelist=seedsAndCores, **seedcorenodeOpts)
		nx.draw_networkx_nodes(G, pos,nodelist=onlyCores, **corenodeOpts)
		nx.draw_networkx_nodes(G, pos, nodelist=others, **nodeOpts)
		
		
		
	def drawEdgeSets(self, G, edgeSets, palette, nodeOpts={}, labelOpts={}, edgeOpts={}, pos=None):
		"""
		Draw different set of edges in different colors

		:param G: The graph object
		:param edgeSets: List of edgesets
		:param palette: List of colors. Every edgeset is drawn in another color of palette.
		:param nodeOpts: Dictionary that contains the node attributes
		:param labelOpts: Dictionary that contains the label attributes
		:param edgeOpts: Dictionary that contains the edge attributes
		"""
		
		plt.figure(figsize=self.size)
		plt.axis("off")
		if not pos:
			pos = self.layout(G)
		
		# draw nodes
		nx.draw_networkx_nodes(G, pos, **nodeOpts)
		# draw labels
		nx.draw_networkx_labels(G, pos, **labelOpts)
		# draw edges
		for i in range(len(edgeSets)):
			nx.draw_networkx_edges(G, pos, edgeSets[i], edge_color=palette[i], **edgeOpts)
		
		
	def drawWeighted(self, G, edge_cmap=plt.cm.Blues, node_color='c', pos=None):
		"""
		Draw a weighted graph so that the more weighted edges appear in a darker color

		:param G: The graph object
		:param edge_cmap: Matplotlib colormap for edgecolors
		:param node_color: The color of nodes
		:param pos: optional predefined node positions
		"""
		
		plt.figure(figsize=self.size)
		plt.axis("off")
		if not pos:
			pos = self.layout(G)

		
		# check if G is weighted
		# if there is at least one weighted edge, the graph is weighted
		weight = False

		for u,v in G.edges():
			if 'weight' in G[u][v]:
				weight = True

	
		# create dictionary for labelling the weighted edges
		# if weight = 1 the edge is unlabelled
		if weight:
			label = {(u,v): G[u][v] for (u,v) in G.edges()}
			if (label[(u,v)] == 1 for (u,v) in G.edges()):
				label = {(u,v): ' ' for (u,v) in G.edges()}

		# create list for colormap		
		colors = []
		if weight:
			colors += [G[u][v]['weight'] for (u,v) in G.edges()]
		if not weight:
			colors += [1 for (u,v) in G.edges()]

		# extra edgeOpts needed for visualising weighted edges	
		edgeOpts = {"edge_color": colors, "edge_cmap": edge_cmap, "edge_vmin": 0, "edge_vmax": 2.5}
		edge_labelOpts = {"alpha": 0.8, "edge_labels": label, "font_size": 10}

		# draw edges
		# TODO: drawWeighted should also work if emphasizeCliques is true
		if self.emphCliques:	
##			alpha = self.emphasizeCliques(G)
##			i = 0
##			for u,v in G.edges():
##			        nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], edge_color=colors,
##						edge_cmap= edge_cmap, edge_vmin=edge_vmin, edge_vmax=edge_vmax,
##						alpha=alpha[i])
##				#i += 1
			pass
		else:
			nx.draw_networkx_edges(G, pos, **edgeOpts)

		# draw nodes
		nx.draw_networkx_nodes(G, pos, **self.nodeOpts)
		
		# draw edgelabels
		if weight:
			nx.draw_networkx_edge_labels(G, pos, **edge_labelOpts)
			
		# draw nodelabels
		if self.labelled:
			nx.draw_networkx_labels(G, pos, **self.labelOpts)
	

		
		
	def drawClustered(self, G, clustering, cmap=plt.cm.Accent, pos=None):
		"""
		Draw a graph with a clustering so that nodes in different clusterings appear in different colors

		:param G: The graph object
		:param pos: optional predefined node positions
		:param clustering: Cluster of graph G
		:param cmap: colormap for visualizing cluster. Default is Accent from Matplotlib
		"""
		if not pos:
			pos = self.layout(G)
		plt.figure(figsize=self.size)
		plt.axis("off")
		
		colors = []
		colors += [clustering[i] for i in clustering]
		
		
		edgeOpts = {}
		nodeOpts = {"node_color": colors, "cmap": cmap} # extra node opts needed for community visualization
		nodeClusterOpts = dict(self.nodeOpts, **nodeOpts)

		# draw nodes
		nx.draw_networkx_nodes(G, pos, **nodeClusterOpts)

		# draw edges
		if self.emphCliques:
			alpha = self.emphasizeCliques(G)
			i = 0
			for u,v in G.edges():
				nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], alpha=alpha[i])
				i += 1
		else:
			nx.draw_networkx_edges(G, pos, **self.edgeOpts)

		# draw nodelabels	
		if self.labelled:
			nx.draw_networkx_labels(G, pos, **self.labelOpts)
		



