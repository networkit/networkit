from warnings import warn

from .graphio import GraphMLSAX as IoGraphMLSAX
from .graphio import GraphMLReader as IoGraphMLReader
from .graphio import GraphMLWriter as IoGraphMLWriter

# GraphML Reader
class GraphMLSAX(IoGraphMLSAX):
	""" 
	GraphMLSAX()

	DEPRECATED. This class (and the networkit.GraphMLIO module) will be removed in future updates.
	Use networkit.graphio.GraphMLSAX instead.

	Parser for GraphML XML files, based on Pythons XML.SAX implementation.
	"""

	def __init__(self):
		""" Initializes several important variables """
		warn("networkit.GraphMLIO.GraphMLSAX is deprecated, will be removed in future updates. Use networkit.graphio.GraphMLSAX instead.")
		super().__init__()

class GraphMLReader(IoGraphMLReader):
	""" 
	GraphMLReader()

	DEPRECATED. This class (and the networkit.GraphMLIO module) will be removed in future updates.
	Use networkit.graphio.GraphMLReader instead.

	This class serves as wrapper for the GraphMLSAX class
	which is able to parse a GraphML XML file and construct
	a graph.
	"""

	def __init__(self):
		warn("networkit.GraphMLIO.GraphMLReader is deprecated, will be removed in future updates. Use networkit.graphio.GraphMLReader instead.")
		super().__init__()

# GraphMLWriter
class GraphMLWriter(IoGraphMLWriter):
	""" 
	GraphMLWriter()

	DEPRECATED. This class (and the networkit.GraphMLIO module) will be removed in future updates.
	Use networkit.graphio.GraphMLWriter instead.

	This class provides a function to write a NetworKit graph to a file in the 
	GraphML format.
	"""
	
	def __init__(self):
		""" Initializes the class. """
		warn("networkit.GraphMLIO.GraphMLWriter is deprecated, will be removed in future updates. Use networkit.graphio.GraphMLWriter instead.")
		super().__init__()
