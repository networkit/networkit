import queue
import xml.etree.cElementTree as ET
from xml.dom import minidom
from _NetworKit import Graph, GraphEvent

# GEXF Reader
class GEXFReader:
	def __init__(self):
		""" Initializes the GEXFReader class """
		self.mapping = dict()
		self.g = Graph(0)
		self.weighted = False
		self.directed = False
		self.dynamic = False
		self.hasDynamicWeights = False
		self.q = queue.Queue()
		self.eventStream = []
		self.nInitialNodes = 0
		self.timeFormat = ""

	def read(self, fpath):
		""" Reads and returns the graph object defined in fpath """
		#0. Reset internal vars and parse the xml
		self.__init__()
		doc = minidom.parse(fpath)

		#1. Determine if graph is dynamic, directed and has dynamically changing weights
		graph = doc.getElementsByTagName("graph")[0]
		if (graph.getAttribute("defaultedgetype") == "directed"):
			self.directed = True
		if (graph.getAttribute("mode") == "dynamic"):
			self.dynamic = True
		if self.dynamic:
			self.timeFormat = graph.getAttribute("timeformat")
		attributes = graph.getElementsByTagName("attribute")
		for att in attributes:
			if att.getAttribute("id") == "weight":
				self.hasDynamicWeights = True
				self.weighted = True

		#2. Read nodes and map them to IDs defined in GEXF file
		nodes = doc.getElementsByTagName("node")
		for n in nodes:
			u = n.getAttribute("id")
			if self.dynamic:
				"""
				A GEXF ID can be a string. However, this version of parser accepts ids
				in only 2 formats:
				1. id = "0,1,2," etc.
				2. id = "n0, n1, n2" etc.
				So either an integer or an integer that has n prefix.
				Gephi generates its random graphs in 2nd format for example.
				"""
				_id = ""
				try:
					_id = int(u)
				except:
					_id = int(u[1:])
				# 2-way mapping to refer nodes back in mapDynamicNodes() method
				self.mapping[u] = _id
				self.mapping[_id] = u
				controlList = {'elementAdded': False, 'elementDeleted': False}
				spells = n.getElementsByTagName("spell")
				if len(spells) > 0:
					for s in spells:
						self.parseDynamics(s, "n", controlList, u)
				else:
					self.parseDynamics(n, "n", controlList, u)
			else:
				self.mapping[u] = self.nInitialNodes
				self.nInitialNodes +=1
		if self.dynamic:
			self.mapDynamicNodes()

		#3. Read edges and determine if graph is weighted
		edges = doc.getElementsByTagName("edge")
		for e in edges:
			u = e.getAttribute("source")
			v = e.getAttribute("target")
			w = "1.0"
			if e.hasAttribute("weight"):
				self.weighted = True
				w = e.getAttribute("weight")
			if self.dynamic:
				controlList = {'elementAdded': False, 'elementDeleted': False}
				spells = e.getElementsByTagName("spell")
				if len(spells) > 0:
					for s in spells:
						self.parseDynamics(s, "e", controlList, u, v, w)
				else:
					self.parseDynamics(e, "e", controlList, u, v, w)
			else:
				self.q.put((u, v, w))

		#4. Create graph object
		self.g = Graph(self.nInitialNodes, self.weighted, self.directed)

		#5. Add initial edges to the graph and sort the eventStream by time
		#5.1 Adding initial edges
		while not self.q.empty():
			edge = self.q.get()
			(u, v, w) = (edge[0], edge[1], float(edge[2]))
			self.g.addEdge(self.mapping[u], self.mapping[v], w)

		#5.2 Sorting the eventStream by time and adding timeStep between events that happen in different times
		self.eventStream.sort(key=lambda x:x[1])
		for i in range(1, len(self.eventStream)):
			if self.eventStream[i][1] != self.eventStream[i-1][1]:
				self.eventStream.append((GraphEvent(GraphEvent.TIME_STEP, 0, 0, 0), self.eventStream[i-1][1]))
		self.eventStream.sort(key=lambda x:x[1])
		self.eventStream = [event[0] for event in self.eventStream]
		return (self.g, self.eventStream)



	def parseDynamics(self, element, elementType, controlList,  u,  v = "0", w = "0"):
		"""
			Determine the operations as follows:
			1.Element has start and not deleted before: Create add event
			2.Element has start and deleted before: Create restore event
			3.Element has end:Create del event
			4.If an element has end before start(or no start at all), add it to the initial graph
			5.For dynamic edges, simply go over the attvalues and create
			weight update events

			* A dynamic element must be defined either using only spells
			or inline attributes. These 2 shouldn't be mixed.
			(For example, Gephi will treat them differently. It'll ignore the inline declaration
			if the same element also contains spells)

		"""
		startTime = element.getAttribute("start")
		if startTime == "":
			startTime = element.getAttribute("startopen")
		endTime	= element.getAttribute("end")
		if endTime == "":
			endTime	= element.getAttribute("endopen")
		if self.timeFormat != "date":
			try:
				startTime = float(startTime)
			except:
				pass
			try:
				endTime = float(endTime)
			except:
				pass
				
		if startTime != "" and endTime != "":
			if startTime < endTime and not controlList['elementDeleted']:
				self.createEvent(startTime, "a"+elementType, u, v, w)
				controlList['elementAdded'] = True
			else:
				self.createEvent(startTime, "r"+elementType, u, v, w)
			self.createEvent(endTime, "d"+elementType, u, v, w)
			controlList['elementDeleted'] = True

		if startTime != "" and endTime == "":
			if controlList['elementDeleted']:
				self.createEvent(startTime, "r"+elementType, u, v, w)
			else:
				self.createEvent(startTime, "a"+elementType, u, v, w)
				controlList['elementAdded'] = True

	 	# Handle dynamic edge weights here
		if elementType == "e" and self.hasDynamicWeights:
			attvalues = element.getElementsByTagName("attvalue")
			# If a spell is traversed, attvalues are siblings
			if len(attvalues) == 0:
				attvalues = element.parentNode.parentNode.getElementsByTagName("attvalue")
			for att in attvalues:
				if att.getAttribute("for") == "weight":
					w = att.getAttribute("value")
					startTime = att.getAttribute("start")
					if startTime == "":
						startTime = att.getAttribute("startopen")
					if self.timeFormat != "date":
						startTime = float(startTime)
					# If this edge is not added, first weight update indicates edge addition
					if not controlList['elementAdded']:
						self.createEvent(startTime, "a"+elementType, u, v, w)
						controlList['elementAdded'] = True
					else:
						self.createEvent(startTime, "c"+elementType, u, v, w)

		if startTime == "":
			if not controlList['elementAdded']:
				if elementType == "n":
					self.mapping[u] = self.nInitialNodes
					self.nInitialNodes += 1
				else:
					self.q.put((u,v,w))
				controlList['elementAdded'] = True
			if endTime != "":
				self.createEvent(endTime, "d"+elementType, u, v, w)
				controlList['elementDeleted'] = True

	def createEvent(self, eventTime, eventType, u, v, w):
		"""
			 Creates a NetworKit::GraphEvent from the supplied parameters
			 and passes it to eventStream
		"""
		event, u = None, self.mapping[u]
		if eventType[1] == "e":
			v, w = self.mapping[v], float(w)
		if eventType == "an":
			event = GraphEvent(GraphEvent.NODE_ADDITION, u, 0, 0)
		elif eventType == "dn":
			event = GraphEvent(GraphEvent.NODE_REMOVAL, u, 0, 0)
		elif eventType == "rn":
			event = GraphEvent(GraphEvent.NODE_RESTORATION, u, 0, 0)
		elif eventType == "ae" or eventType == "re":
			event = GraphEvent(GraphEvent.EDGE_ADDITION, u, v, w)
		elif eventType == "de":
			event = GraphEvent(GraphEvent.EDGE_REMOVAL, u, v, w)
		elif eventType == "ce":
			event = GraphEvent(GraphEvent.EDGE_WEIGHT_UPDATE, u, v, w)
		self.eventStream.append((event, eventTime))

	def mapDynamicNodes(self):
		"""
			Node ID of a dynamic node must be determined before it's mapped to its GEXF ID.
			This requires processing the sorted eventStream and figuring out the addition order of the nodes.
			After that, node addition/deletion/restoration operations of this node must be readded to eventStream
			with correct mapping.

			!Note: New mapping of a node can be equal to old mapping of a node. In order to prevent collisions,
			isMapped array must be maintained and controlled.
		"""
		nNodes = self.nInitialNodes
		nEvent = len(self.eventStream)
		isMapped = [False] * nEvent
		self.eventStream.sort(key=lambda x:x[1])
		for i in range(0, nEvent):
			event = self.eventStream[i]
			# Only the nodes with addition event will get remapped.
			if not isMapped[i] and event[0].type == GraphEvent.NODE_ADDITION:
				u = event[0].u
				self.mapping[self.mapping[u]] = nNodes
				# All the other events of that node comes after it's addition event
				for j in range(i, len(self.eventStream)):
					event = self.eventStream[j]
					if not isMapped[j] and event[0].u == u:
						mappedEvent = GraphEvent(event[0].type, self.mapping[self.mapping[u]], 0, 0)
						self.eventStream[j] = (mappedEvent, event[1])
						isMapped[j] = True
				nNodes +=1
				isMapped[i] = True

	def getNodeMap(self):
		""" Returns GEXF ID -> NetworKit::Graph node ID mapping. """
		forwardMap = dict()
		for key in self.mapping:
			if type(key) == str:
				forwardMap[key] = self.mapping[key]
		return forwardMap


# GEXFWriter
class GEXFWriter:
	""" This class provides a function to write a NetworKit graph to a file in the
		GEXF format. """

	def __init__(self):
		""" Initializes the class. """
		self.edgeIdctr = 0
		self.q = queue.Queue()
		self.hasDynamicWeight = False

	def write(self, graph, fname, eventStream = [], mapping = []):
		"""
			Writes a NetworKit::Graph to the specified file fname.
			Parameters:
				- graph: a NetworKit::Graph python object
				- fname: the desired file path and name to be written to
				- eventStream: stream of events
				- mapping: random node mapping
		"""
		#0. Reset internal vars
		self.__init__()

		#1. Start with the root element and the right header information
		root = ET.Element('gexf')
		root.set("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
		root.set("xsi:schemaLocation","http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd")
		root.set('version', '1.2')

		#2. Create graph element with appropriate information
		graphElement = ET.SubElement(root,"graph")
		if graph.isDirected():
			graphElement.set('defaultedgetype', 'directed')
		else:
			graphElement.set('defaultedgetype', 'undirected')
		if len(eventStream) > 0:
			graphElement.set('mode', 'dynamic')
			graphElement.set('timeformat', 'double')
			for event in eventStream:
				if event.type == GraphEvent.EDGE_WEIGHT_UPDATE:
					dynamicAtt = ET.SubElement(graphElement, "attributes")
					dynamicAtt.set('class', 'edge')
					dynamicAtt.set('mode', 'dynamic')
					dynamicWeight = ET.SubElement(dynamicAtt, "attribute")
					dynamicWeight.set('id', 'weight')
					dynamicWeight.set('title', 'Weight')
					dynamicWeight.set('type', 'float')
					self.hasDynamicWeight = True
					break
		else:
			graphElement.set('mode', 'static')

		#3. Add nodes
		nodesElement = ET.SubElement(graphElement, "nodes")
		nNodes, idArray = 0, []
		#3.1 Count the # of nodes (inital + dynamic nodes)
		for event in eventStream:
			if event.type == GraphEvent.NODE_ADDITION:
				nNodes +=1
		nNodes += len(graph.nodes())
		for i in range(0, nNodes):
			idArray.append(i)
		# Optional:Map nodes to a random mapping if user provided one
		if (len(mapping) > 0):
			if(nNodes != len(mapping)):
				 raise Exception('Size of nodes and mapping must match')
			else:
				for i in range(0, nNodes):
					idArray[i] = mapping[i]

		#3.2 Write nodes to the gexf file
		for n in range(nNodes):
			nodeElement = ET.SubElement(nodesElement,'node')
			nodeElement.set('id', str(idArray[n]))
			self.writeEvent(nodeElement, eventStream, n)

		#4. Add edges
		edgesElement = ET.SubElement(graphElement, "edges")
		#4.1 Put all edges into a queue(inital + dynamic edges)
		for e in graph.edges():
			self.q.put((e[0], e[1], graph.weight(e[0], e[1])))
		for event in eventStream:
			if event.type == GraphEvent.EDGE_ADDITION:
				self.q.put((event.u, event.v, event.w))
		#4.2 Write edges to the gexf file
		while not self.q.empty():
			edgeElement = ET.SubElement(edgesElement,'edge')
			e = self.q.get()
			edgeElement.set('source', str(idArray[e[0]]))
			edgeElement.set('target', str(idArray[e[1]]))
			edgeElement.set('id', "{0}".format(self.edgeIdctr))
			self.edgeIdctr += 1
			if graph.isWeighted():
				edgeElement.set('weight', str(e[2]))
			self.writeEvent(edgeElement, eventStream, e)

		#5. Write the generated tree to the file
		tree = ET.ElementTree(root)
		tree.write(fname,"utf-8",True)

	def writeEvent(self, xmlElement, eventStream, graphElement):
		# A var that indicates if the event belongs the graph element we traverse on
		matched = False
		startEvents = [GraphEvent.NODE_ADDITION, GraphEvent.EDGE_ADDITION, GraphEvent.NODE_RESTORATION]
		endEvents = [GraphEvent.NODE_REMOVAL, GraphEvent.EDGE_REMOVAL]
		nodeEvents = [GraphEvent.NODE_ADDITION, GraphEvent.NODE_REMOVAL, GraphEvent.NODE_RESTORATION]
		edgeEvents = [GraphEvent.EDGE_ADDITION, GraphEvent.EDGE_REMOVAL, GraphEvent.EDGE_WEIGHT_UPDATE]
		spellTag, weightTag, operation = False, False, ""
		timeStep = 0
		spellsElement, attValuesElement = None, None

		for event in eventStream:
			if event.type == GraphEvent.TIME_STEP:
				timeStep += 1
			if type(graphElement) == type(0): #a node is an integer
				matched = (event.type in nodeEvents and event.u == graphElement)
			else:
				matched = (event.type in edgeEvents and (event.u == graphElement[0] and event.v == graphElement[1]))
			if matched:
				# Handle weight update seperately
				if event.type == GraphEvent.EDGE_WEIGHT_UPDATE:
					if not weightTag:
						attvaluesElement = ET.SubElement(xmlElement, "attvalues")
						weightTag = True
					attvalue = ET.SubElement(attvaluesElement, "attvalue")
					attvalue.set('for', 'weight')
					attvalue.set('value', str(event.w))
					attvalue.set('start', str(timeStep))
					attvalue.set('endopen', str(timeStep + 1))
				else:
					if event.type in startEvents:
						operation = "start"
					else:
						operation = "end"
					if not spellTag:
						spellsElement = ET.SubElement(xmlElement, "spells")
						spellTag = True
					spellElement = ET.SubElement(spellsElement, "spell")
					spellElement.set(operation, str(timeStep))
