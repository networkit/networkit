from warnings import warn

from .graphio import GEXFReader as IoGEXFReader
from .graphio import GEXFWriter as IoGEXFWriter

# GEXF Reader
class GEXFReader(IoGEXFReader):
	"""
	GEXFReader()

	DEPRECATED. This class (and the networkit.GEXFIO module) will be removed in future updates.
	Use networkit.graphio.GEXFReader instead.

	This class provides a function to read a file in the
	GEXF (Graph Exchange XML Format) format.

	For more details see: http://gexf.net/
	"""
	def __init__(self):
		""" Initializes the GEXFReader class """
		warn("networkit.GEXFIO.GEXFReader is deprecated, will be removed in future updates. Use networkit.graphio.GEXFReader instead.")
		super().__init__()

# GEXFWriter
class GEXFWriter(IoGEXFWriter):
	""" 
	GEXFWriter()

	DEPRECATED. This class (and the networkit.GEXFIO module) will be removed in future updates.
	Use networkit.graphio.GEXFWriter instead.
	
	This class provides a function to write a NetworKit graph to a file in the GEXF format. 
	"""

	def __init__(self):
		""" Initializes the class. """
		warn("networkit.GEXFIO.GEXFWriter is deprecated, will be removed in future updates. Use networkit.graphio.GEXFWriter instead.")
		super().__init__()
