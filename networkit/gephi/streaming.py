"""
This module provides methods to export data from NetworKit directly to Gephi using the
Gephi Streaming Plugin.
"""

__author__ = "Gerd Lindner, Moritz v.Looz"


import urllib as _urllib
import time
import math

from . import pyclient as _gephipyclient   # we want to hide these from the user

class GephiStreamingClient:
	"""
	GephiStreamingClient(url='http://localhost:8080/workspace0')

	A python singleton providing access to a running Master instance of the Gephi
	Streaming plugin.

	Parameters
	----------
	url : str, optional
		URL of the Gephi server. Default: http://127.0.0.1:8080/workspace0
	"""

	def __init__(self, url='http://localhost:8080/workspace0'):
		#Disabling Flushing means quite a good performance boost.
		self._pygephi = _gephipyclient.GephiClient(url, autoflush=10000)
		self.graphExported = False

	def _urlError(self, e):
		print("Could not connect to the gephi streaming plugin.")
		print("Did you start the streaming master server in gephi and provide the name of your workspace?")
		print("If the workspace is named 'Workspace 0', the corresponding url is http://localhost:8080/workspace0 (adapt port)")

	def _edgeId(self, u, v):
		""" Return the edge id that is exported to Gephi
		"""
		if self.directed:
			return str(u) + '->' + str(v)
		return str(min(u,v)) + '-' + str(max(u,v))

	def exportGraph(self, graph):
		""" 
		exportGraph(graph)
		
		Exports the given graph to gephi. No attributes or weights are exported.
		Please note that only graphs with indexed edges can be exported, so
		edges will be indexed if neccessary. 
		
		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		"""

		try:
			self._pygephi.clean()

			#Index edges if neccessary
			graph.indexEdges()
			self.directed = graph.isDirected()

			self._exportNodes(graph)

			for u, v in graph.iterEdges():
				self._pygephi.add_edge(self._edgeId(u, v), u, v, self.directed)

			self._pygephi.flush()
			self.graphExported = True
		except _urllib.error.URLError as e:
			self._urlError(e)

	def _exportNodes(self, graph):
		nAttrs = {'size': 2.0, 'r': 0.6, 'g': 0.6, 'b': 0.6, 'y':1.0}

		# the default approximately shows -2000 to 2000, so we want to
		# distribute the nodes in that area. Since Gephi 0.9, no nodes
		# may have exactly the same coordinates, thus a deterministic
		# distribution scheme is used.
		NODE_AREA_SIZE = 2000
		nodesPerSquareSide = 0 if graph.numberOfNodes() == 0 else math.ceil(math.sqrt(graph.numberOfNodes()))
		stepSize = NODE_AREA_SIZE / nodesPerSquareSide
		offset = NODE_AREA_SIZE / 2

		for nodeNumber, node in enumerate(graph.iterNodes()):
			nAttrs['x'] = (nodeNumber % nodesPerSquareSide) * stepSize - offset
			nAttrs['y'] = (nodeNumber // nodesPerSquareSide) * stepSize - offset
			self._pygephi.add_node(str(node), **nAttrs)

	def exportAdditionalEdge(self, u, v):
		"""
		exportAdditionalEdge(u, v)
		
		Adds an edge (u,v) in an already exported graph. If the edge is already present, nothing happens.
		If the graph is directed, the edge goes from u to v.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.
		"""
		if self.graphExported != True:
			print("Error: Cannot add edges. Export Graph first!")
			return
		try:
			self._pygephi.add_edge(self._edgeId(u, v), u, v, self.directed)
			self._pygephi.flush()
		except _urllib.error.URLError as e:
			self._urlError(e)

	def removeExportedEdge(self, u, v):
		""" 
		removeExportedEdge(u, v)
		
		Removes an edge from an already exported graph.
		
		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.
		"""
		if self.graphExported != True:
			print("Error: Cannot remove edges. Export Graph first!")
			return
		try:
			self._pygephi.delete_edge(self._edgeId(u, v))
			self._pygephi.flush()
		except _urllib.error.URLError as e:
			self._urlError(e)

	def exportEventStream(self, stream, timeStepDelay = 0):
		"""
		exportEventStream(stream, timeStepDelay = 0)

		Export a given event stream of networkit.dynamics.GraphEvent to Gephi server.

		Parameters
		----------
		stream : list(networkit.dynamics.GraphEvent)
			The input stream.
		timeStepDelay : int, optional
			Add a delay (in ms) for every time step event. Default: 0
		"""
		if self.graphExported != True:
			print("Error: Cannot export event stream. Export Graph first!")
			return

		nAttrs = {}
		try:
			for ev in stream:
				if ev.type == ev.NODE_ADDITION:
					self._pygephi.add_node(str(ev.u), **nAttrs)
				elif ev.type == ev.NODE_REMOVAL:
					self._pygephi.delete_node(str(ev.u))
				elif ev.type == ev.EDGE_ADDITION:
					self._pygephi.add_edge(self._edgeId(ev.u, ev.v), ev.u, ev.v, self.directed)
				elif ev.type == ev.EDGE_REMOVAL:
					self._pygephi.delete_edge(self._edgeId(ev.u, ev.v))
				elif ev.type == ev.EDGE_WEIGHT_UPDATE:
					print("Edge weights not yet supported in gephi streaming!")
				elif ev.type == ev.EDGE_WEIGHT_INCREMENT:
					print("Edge weights not yet supported in gephi streaming!")
				elif ev.type == ev.TIME_STEP:
					self._pygephi.flush()
					if timeStepDelay > 0:
						time.sleep(timeStepDelay)
				self._pygephi.flush()
		except _urllib.error.URLError as e:
			self._urlError(e)

	def exportNodeValues(self, graph, values, attribute_name):
		"""
		exportNodeValues(graph, values, attribute_name)

		This method exports node values (e.g. community information, betwenness centrality values)
		to gephi using the Gephi Streaming Plugin. Use exportGraph(Graph) first to export
		the graph itself.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		values : list(float) or networkit.Partition
			Python list or Partition object that contains the values to be exported.
		attribute_name : str
			Name of the node attribute in Gephi.
		"""

		try:
			if len(values) != graph.numberOfNodes():
				print("Warning: #Nodes (", graph.numberOfNodes(), ") does not match #Values (", len(values), ").")

			for i in graph.iterNodes():
				nAttrs = {attribute_name:values[i]}
				self._pygephi.change_node(str(i), **nAttrs)

			self._pygephi.flush()
		except _urllib.error.URLError as e:
			self._urlError(e)

	def exportCoordinates(self, graph, scale=1):
		try:
			xcoords = [scale*graph.getCoordinate(v)[0] for v in graph.iterNodes()]
			ycoords = [scale*graph.getCoordinate(v)[1] for v in graph.iterNodes()]
			self.exportNodeValues(graph, xcoords, 'x')
			self.exportNodeValues(graph, ycoords, 'y')
			self._pygephi.flush()
		except _urllib.error.URLError as e:
			self._urlError(e)

	def exportEdgeValues(self, graph, values, attribute_name):
		"""
		exportEdgeValues(graph, values, attribute_name)

		This method exports an edge attribute to gephi using the Gephi Streaming Plugin.
		Use exportGraph(Graph) first to export graph itself.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		values : list(float) or networkit.Partition
			Python list or Partition object that contains the values to be exported.
		attribute_name : str
			Name of the node attribute in Gephi.
		"""
		try:
			if len(values) != graph.upperEdgeIdBound():
				print("Warning: Upper edge id bound (", graph.upperEdgeIdBound(), ") does not match #Values (", len(values), ").")

			edgetype = "Directed" if self.directed else "Undirected"

			def exportEdgeV(u, v, _w, eid):
				eAttrs = {attribute_name:values[eid], "Type":edgetype}
				self._pygephi.change_edge(self._edgeId(u, v), u, v, self.directed, **eAttrs)

			graph.forEdges(exportEdgeV)

			self._pygephi.flush()
		except _urllib.error.URLError as e:
			self._urlError(e)

	def clearGraph(self):
		""" 
		clearGraph()
		
		Cleans the first workspace in Gephi.
		"""

		try:
			self._pygephi.clean()
			self._pygephi.flush()
			self.graphExported = False

		except _urllib.error.URLError as e:
			self._urlError(e)
