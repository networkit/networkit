import xml.etree.cElementTree as ET
import xml.sax

from _NetworKit import Graph

# GraphML Reader
class GraphMLSAX(xml.sax.ContentHandler):
	""" Parser for GraphML XML files, based on Pythons XML.SAX implementation. """

	def __init__(self):
		""" Initializes several important variables """
		xml.sax.ContentHandler.__init__(self)
		self.charBuffer = []
		self.mapping = dict()
		self.g = Graph(0)
		self.graphName = 'unnamed'
		self.weightedID = ''
		self.weighted = False
		self.directed = False
		self.edgestack = []
		self.edgeweight = 0.0
		self.keepData = False

	def startElement(self, name, attrs):
		""" Parses all currently relevant XML tags and retrieves data."""
		if name == "graph":
			# determine, if graph is directed:
			if attrs.getValue("edgedefault") == "directed":
				print("identified graph as directed")
				self.directed = True
			if "id" in  attrs.getNames() and not attrs.getValue("id") == '':
					self.graphName = attrs.getValue("id")
			self.g = Graph(0,self.weighted, self.directed)
			self.g.setName(self.graphName)
		if name == "node":
			u = self.g.addNode()
			val = attrs.getValue("id")
			self.mapping[val] = u
		elif name == "edge":
			u = attrs.getValue("source")
			v = attrs.getValue("target")
			self.edgestack.append((u,v))
		elif name == "key":
			#print("found element with tag KEY")
			if (attrs.getValue("for") == 'edge' and attrs.getValue("attr.name") == 'weight' and attrs.getValue("attr.type") == 'double'):
				self.weighted = True
				self.weightedID = attrs.getValue("id")
				print("identified graph as weighted")
		elif name == "data" and attrs.getValue("key") == self.weightedID:
			self.keepData = True

	def endElement(self, name):
		""" Finalizes parsing of the started Element and processes retrieved data."""
		data = self.getCharacterData()
		if name == "edge":
			u = self.edgestack[len(self.edgestack)-1][0]
			v = self.edgestack[len(self.edgestack)-1][1]
			self.edgestack.pop()
			if self.weighted:
				#print ("identified edge as weighted with weight: {0}".format(edgeweight))
				self.g.addEdge(self.mapping[u], self.mapping[v], self.edgeweight)
				self.edgeweight = 0.0
			else:
				self.g.addEdge(self.mapping[u], self.mapping[v])
		elif name == "data" and self.keepData:
			self.keepData = False
			self.edgeweight = float(data)

	def characters(self, content):
		self.charBuffer.append(content)

	def getCharacterData(self):
		data = ''.join(self.charBuffer).strip()
		self.charBuffer = []
		return data

	def getGraph(self):
		return self.g

class GraphMLReader:
	""" This class serves as wrapper for the GraphMLSAX class
		which is able to parse a GraphML XML file and construct
		a graph. """

	def __init__(self):
		""" Initializes the GraphMLSAX class """
		self.graphmlsax = GraphMLSAX()

	def read(self, fpath):
		""" Parses a GraphML XML file and returns the constructed Graph
			Parameters:
				- fpath: the path to the file as a string
		"""
		xml.sax.parse(fpath, self.graphmlsax)
		return self.graphmlsax.getGraph()

# GraphMLWriter
class GraphMLWriter:
	""" This class provides a function to write a NetworKit graph to a file in the 
		GraphML format. """
	
	def __init__(self):
		""" Initializes the class. """
		self.edgeIdCounter = 0
		self.dir_str = ''

	def write(self, graph, fname, nodeAttributes = {}, edgeAttributes = {}):
		""" Writes a NetworKit graph to the specified file fname. 
			Parameters:
				- graph: a NetworKit::Graph python object 
				- fname: the desired file path and name to be written to
				- nodeAttributes: optional dictionary of node attributes in the form attribute name => list of attribute values
				- edgeAttributes: optional dictionary of edge attributes in the form attribute name => list of attribute values
		"""
		# reset some internal variables in case more graphs are written with the same instance
		self.edgeIdCounter = 0
		self.dir_str = ''

		if len(edgeAttributes) > 0 and not graph.hasEdgeIds():
			raise RuntimeError("Error, graph must have edge ids if edge attributes are given")

		# start with the root element and the right header information
		root = ET.Element('graphml')
		root.set("xmlnsi","http://graphml.graphdrawing.org/xmlns")
		root.set("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
		root.set("xsi:schemaLocation","http://graphml.graphdrawing.org/xmlns \
			http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd")

		maxAttrKey = 1
		# if the graph is weighted, add the attribute
		if graph.isWeighted():
			attrElement = ET.SubElement(root,'key')
			attrElement.set('for','edge')
			attrElement.set('id', 'd1')
			attrElement.set('attr.name','weight')
			attrElement.set('attr.type','double')
			maxAttrKey += 1

		attrKeys = {}
		import numbers
		import itertools
		for attType, attName, attData in itertools.chain(
			map(lambda d : ('node', d[0], d[1]), nodeAttributes.items()),
			map(lambda d : ('edge', d[0], d[1]), edgeAttributes.items())):

			attrElement = ET.SubElement(root, 'key')
			attrElement.set('for', attType)
			attrElement.set('id', 'd{0}'.format(maxAttrKey))
			attrKeys[(attType, attName)] = 'd{0}'.format(maxAttrKey)
			maxAttrKey += 1
			attrElement.set('attr.name', attName)
			if len(attData) == 0:
				attrElement.set('attr.type', 'int')
			elif isinstance(attData[0], bool):
				attrElement.set('attr.type', 'boolean')
				# special handling for boolean attributes: convert boolean into lowercase string
				if attType == 'edge':
					edgeAttributes[attName] = [ str(d).lower() for d in attData ]
				else:
					nodeAttributes[attName] = [ str(d).lower() for d in attData ]
			elif isinstance(attData[0], numbers.Integral):
				attrElement.set('attr.type', 'int')
			elif isinstance(attData[0], numbers.Real):
				attrElement.set('attr.type', 'double')
			else:
				attrElement.set('attr.type', 'string')


		# create graph element with appropriate information
		graphElement = ET.SubElement(root,"graph")
		if graph.isDirected():
			graphElement.set('edgedefault', 'directed')
			self.dir_str = 'true'
		else:
			graphElement.set('edgedefault', 'undirected')
			self.dir_str = 'false'
		graphElement.set('id', graph.getName())

		# Add nodes
		for n in graph.nodes():
			nodeElement = ET.SubElement(graphElement,'node')
			nodeElement.set('id', str(n))
			for attName, attData in nodeAttributes.items():
				dataElement = ET.SubElement(nodeElement, 'data')
				dataElement.set('key', attrKeys[('node', attName)])
				dataElement.text = str(attData[n])

		# in the future: more attributes
	        #for a in n.attributes():
        	#    if a != 'label':
	        #        data = doc.createElement('data')
        	#        data.setAttribute('key', a)
	        #        data.appendChild(doc.createTextNode(str(n[a])))
        	#        node.appendChild(data)

		# Add edges
		def addEdge(u, v, w, eid):
			edgeElement = ET.SubElement(graphElement,'edge')
			edgeElement.set('directed', self.dir_str)
			edgeElement.set('target', str(v))
			edgeElement.set('source', str(u))
			if graph.hasEdgeIds():
				edgeElement.set('id', "e{0}".format(eid))
			else:
				edgeElement.set('id', "e{0}".format(self.edgeIdCounter))
				self.edgeIdCounter += 1
			if graph.isWeighted():
				# add edge weight
				dataElement = ET.SubElement(edgeElement,'data')
				dataElement.set('key','d1')
				dataElement.text = str(w)
			for attName, attData in edgeAttributes.items():
				dataElement = ET.SubElement(edgeElement, 'data')
				dataElement.set('key', attrKeys[('edge', attName)])
				dataElement.text = str(attData[eid])
		graph.forEdges(addEdge)

	#TODO: optional prettify function for formatted output of xml files
		tree = ET.ElementTree(root)
		tree.write(fname,"utf-8",True)


