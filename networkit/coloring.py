# local imports
from . import graph
#import .
from .algebraic import adjacencyEigenvectors

class SpectralColoring(object):
	"""
	SpectralColoring(G)

	Algorithm for computing a valid spectral coloring for a given graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to compute the spectral coloring for.
	"""
	def __init__(self, G):
		super(SpectralColoring, self).__init__()

		self.graph = G

	def prepareSpectrum(self):
		"""
		prepareSpectrum()

		Computes eigenvalues and eigenvector in order to prepare the spectral coloring. 
		This function does not need to be called, it is sufficient to call run()
		for the SpectralColoring-object.
		"""
		spectrum = adjacencyEigenvectors(self.graph)
		self.eigenvalues = spectrum[0]
		self.eigenvectors = spectrum[1]

	def valid(self, color):
		"""
		valid(color)

		Checks whether colorID is valid

		Parameters
		----------
		color : int
			colorID to check for.

		Returns
		-------
		bool
			Returns True if valid, False if not.
		"""
		for v in self.colors[color]:
			for u in self.graph.iterNeighbors(v):
				if u in self.colors[color]:
					return False

		return True

	def split(self, color, depth=0):
		"""
		split(color, depth=0)

		Helper method for computing spectral coloring. 
		This function does not need to be called, it is sufficient to call run()
		for the SpectralColoring-object.

		Parameters
		----------
		color : int
			colorID to split corresponding nodes into colorID and new color.
		depth : int, optional
			Sets index of eigenvector to compare to `depth`. Default: 0
		"""
		otherColor = self.nextColor
		self.nextColor += 1

		vs = self.colors[color]

		self.colors[color] = [v for v in vs if self.eigenvectors[depth][v] >= 0]
		self.colors[otherColor] = [v for v in vs if self.eigenvectors[depth][v] < 0]

		if not self.valid(color):
			self.split(color, depth=depth+1)

		if not self.valid(otherColor):
			self.split(otherColor, depth=depth+1)

	def buildReverseDict(self):
		"""
		buildReverseDict()

		Helper method for computing spectral coloring. 
		This function does not need to be called, it is sufficient to call run()
		for the SpectralColoring-object.
		"""
		self.coloring = {}

		for color in self.colors:
			for v in self.colors[color]:
				self.coloring[v] = color


	def run(self):
		"""
		run()

		Main method of SpectralColoring. This computes a valid coloring.
		"""
		self.prepareSpectrum()

		self.colors = {0 : set(self.graph.iterNodes())}
		self.nextColor = 1

		self.split(0)
		self.buildReverseDict()

	def getColoring(self):
		"""
		getColoring()

		Returns a valid coloring for a given graph. Each node has a color ID to be mapped onto a color palette.

		Returns
		-------
		list(int)
			A list of nodes containing a valid color mapping. 
		"""
		return self.coloring
