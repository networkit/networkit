# local imports
from . import graph
#import .
from .algebraic import adjacencyEigenvectors

class SpectralColoring(object):
	def __init__(self, G):
		super(SpectralColoring, self).__init__()

		self.graph = G

	def prepareSpectrum(self):
		spectrum = adjacencyEigenvectors(self.graph)
		self.eigenvalues = spectrum[0]
		self.eigenvectors = spectrum[1]

	def valid(self, color):
		for v in self.colors[color]:
			for u in self.graph.neighbors(v):
				if u in self.colors[color]:
					return False

		return True

	def split(self, color, depth=0):
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
		self.coloring = {}

		for color in self.colors:
			for v in self.colors[color]:
				self.coloring[v] = color


	def run(self):
		self.prepareSpectrum()

		self.colors = {0 : set(self.graph.nodes())}
		self.nextColor = 1

		self.split(0)
		self.buildReverseDict()

	def getColoring(self):
		return self.coloring
