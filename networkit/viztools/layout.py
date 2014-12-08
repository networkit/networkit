'''
Created on 02.05.2013

@author: cls
'''
# local imports
from ..algebraic import laplacianEigenvectors

# external imports
import random
import numpy as np
import itertools
import math



class Layout:
	"""Superclass for all layout algorithms"""
	
	def layout(self, G):
		raise NotImplementedError("abstract method")
	
class SpectralLayout(Layout):

	def _prepareSpectrum(self):
		spectrum = laplacianEigenvectors(self.graph, cutoff = 2, reverse=True)
		self.eigenvectors = spectrum[1]
		self.eigenvalues = spectrum[0]

	def _prepareNormalization(self):
		x = self.eigenvectors[1] / self.eigenvalues[1]
		y = self.eigenvectors[2] / self.eigenvalues[2]

		minX = x[0]
		maxX = x[0]

		for val in x:
			minX = min(minX, val)
			maxX = max(maxX, val)

		self.minX = minX
		self.rangeX = maxX - minX

		minY = y[0]
		maxY = y[0]

		for val in y:
			minY = min(minY, val)
			maxY = max(maxY, val)

		self.minY = minY
		self.rangeY = maxY - minY

	def _computePositions(self):
		self.pos = {}

		for v in self.graph.nodes():
			x = ((self.eigenvectors[1][v] / self.eigenvalues[1]) - self.minX) / self.rangeX
			y = ((self.eigenvectors[2][v] / self.eigenvalues[2]) - self.minY) / self.rangeY

			self.pos[v] = (x,y)

	def _getPos(self, v):
		return np.array(self.pos[v])

	def layout(self, G):
		"""
		Perform a spectral layout of the Graph G.

		Parameters
		----------
		G : graph
		    The graph to layout

		Returns
		-------
		A dictionary mapping vertices to 2D positions.
		"""
		self.graph = G
		self._prepareSpectrum()
		self._prepareNormalization()
		self._computePositions()

		return {v : self._getPos(v) for v in self.graph.nodes()}

