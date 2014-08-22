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
import networkx as nx



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


class ForceDirected(Layout):

	def __init__(self, repConst=0.4, step=0.05, scale=1, iterations=30):
		"""
		Formula of attractive force: $\frac{distance}{springLength}$
		Formula of repelling force: $\frac{repConst}{distance}$

		:param repConst: constant of repelling force. The higher repConst the higher becomes the repelling force
		:param step: step determines how much the nodes are moved in each iteration
		:param scale: coordinates are in interval [0,scale]$\times$[0, scale]
		:param iterations: number of iterations
		"""
		self.attractors = []
		self.criterion = []
		self.strength = []
		self.coord = None
		self.repConst = repConst
		self.step = step
		self.scale = scale
		self.iterations = iterations

	def addAttractor(self, attractor=[], criterion=[], strength=10):
		"""An attractor is a point on the canvas which attracts (or repels)
		nodes matching a certain criterion.
		
		:param attractor: Coordinates of the attractor
		:param criterion: Function with a node as parameter. If criterion returns True the node is attracted (or repelled).
				  Example: criterion=lambda v: G.degree(v) == 0
		:param strength: Strength of the acting force of the attractor (default = 10).
                                 The force is repelling if strength is negativ.
		"""
		
		# if attractor, criterion or strength =[] layout does not change
		if attractor==[] or criterion==[] or strength==[]:
			pass
			
		else:
			self.attractors += [attractor]
			self.criterion += [criterion]
			self.strength += [strength]
				
	def layout(self, G):
	
		""" Create a layout (node coordinates) for the graph

		:param G: The graph object
		"""

		# coordinates are in interval [0, 1]
		self.coord = {v : np.array((random.random(), random.random())) for v in G.nodes()}

		
		# In Python you can define functions inside of functions, which is very handy sometimes.
		# In these nested functions you can also use variables from the surrounding scope.
		
		
		
		def distance(p1, p2):
			""" Return the distance between two points"""
			(x1, y1) = p1
			(x2, y2) = p2
			distance = math.sqrt((x1-x2)**2+(y1-y2)**2)
			return distance
			
		
		def attractiveForce(p1, p2):
			""" Return the attractive force of p2 acting on p1"""
			d = distance(p1, p2)
			if d < 0.000001: # Error if d is too small
				d = 0.000001
			springLength = (self.scale**2)/(nx.number_of_nodes(G))
			value = d/springLength
			direction = p2 - p1 
			force = (direction/d) * value # = unit vector * value of the force
			return force
	
		
		
		def repellingForce(p1, p2):
			""" Return the repelling force of p2 acting on p1"""
			d = distance(p1, p2)
			if d < 0.000001: # Error if d is too small
				d = 0.000001
			value = self.repConst/d       
			direction = p1 - p2
			force = (direction/d) * value # = unit vector * value of the force
			return force

		def attractorForce(G):
			""" Return the force of the attractors acting on the nodes of the graph """
			forces = {v: np.array((0.0, 0.0)) for v in G.nodes()}
			
			i = 0
			while i < len(self.attractors):
				
				#force is attractive
				if self.strength[i] > 0:
					for u in G.nodes():
						if self.criterion[i](u):
							forces[u] += attractiveForce(self.coord[u], self.attractors[i])*self.strength[i]
						
				#force is repelling		
				elif self.strength[i] < 0:
					for u in G.nodes():
						if self.criterion[i](u):                                                
							forces[u] += repellingForce(self.coord[u], self.attractors[i])*(-1)*self.strength[i]
				i += 1
			
			return forces
	
		def move(p, force):
			""" Move the point along the direction of the force"""
			
			p = p + force*self.step
			return p

		def rescaleLayout(pos):
			""" Rescale to (0,scale) in all axes
			Source: http://networkx.github.io/documentation/networkx-1.8/_modules/networkx/drawing/layout.html
			"""
			
			# shift origin to (0,0)
			lim = 0 # max coordinate for all axes
			for i in range(2):
				for u in range(nx.number_of_nodes(G)):
					pos[u][i] -= min(pos[u][i] for u in range(nx.number_of_nodes(G)))
				lim=max(max(pos[u][i] for u in range(nx.number_of_nodes(G))),lim)
				
			# rescale to (0,scale) in all directions, preserves aspect
			for i in range(2):
				for u in range(nx.number_of_nodes(G)):
					pos[u][i]*=self.scale/lim
			return pos

		M = 0			
		def done():
			""" Return whether layout is done.
			After a number of iterations layout is done."""
			if M > self.iterations:
				return True
			else:
				return False
				
		# main algorithm
		
		# a force vector acts on every node
		forces = {v: np.array((0.0, 0.0)) for v in G.nodes()}

		while not done():

			M += 1
			
			for (u, v) in itertools.combinations(G.nodes(), 2): # all node pairs
                                
				# 1) repelling forces between nodes
				forces[u] += repellingForce(self.coord[u], self.coord[v])
				forces[v] += repellingForce(self.coord[v], self.coord[u])

			for (u, v) in G.edges():
				# 2) attractive forces between connected nodes
				if (u != v): #self-loops do not create force
					forces[u] += attractiveForce(self.coord[u], self.coord[v])
					forces[v] += attractiveForce(self.coord[v], self.coord[u])

				# 3) optional user-defined forces, e.g. forces that separate classes of nodes

			# optional attractors
			if self.attractors:				
				for u in G.nodes():
					forces[u] += attractorForce(G)[u]

			# move node coordinates
			for u in G.nodes():
				self.coord[u] = move(self.coord[u], forces[u])

			# rescale coordinates to [0, scale]
			self.coord =rescaleLayout(self.coord)
			

		return self.coord

			



