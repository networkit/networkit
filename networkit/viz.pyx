# distutils: language=c++

from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from .helpers import stdstring

from .graph cimport _Graph, Graph
from .structures cimport count, index, node, coordinate

cdef extern from "<networkit/viz/Point.hpp>" namespace "NetworKit" nogil:

	cdef cppclass Point[T]:
		Point()
		Point(T x, T y)
		count getDimensions()
		T& operator[](const index i) except +
		T& at(const index i) except +

	cdef cppclass _Point2D "NetworKit::Point2D":
		_Point2D()
		pair[coordinate, coordinate] asPair()

cdef object toPoint2DVector(const vector[_Point2D]& v):
	return [v[i].asPair() for i in range(v.size())]

cdef object toNodePoint2DVector(const vector[pair[node, _Point2D]]& v):
	return [(v[i].first, v[i].second.asPair()) for i in range(v.size())]

cdef extern from "<networkit/viz/GraphLayoutAlgorithm.hpp>":

	cdef cppclass _GraphLayoutAlgorithm "NetworKit::GraphLayoutAlgorithm"[T]:
		_GraphLayoutAlgorithm(_Graph, count) except +
		count numEdgeCrossings() except +
		vector[Point[double]] getCoordinates() except +
		bool_t writeGraphToGML(string path) except +
		bool_t writeKinemage(string path) except +
		void run() nogil except +

cdef class GraphLayoutAlgorithm:
	"""
	Abstract base class for graph drawing algorithms.
	"""

	cdef _GraphLayoutAlgorithm[double] *_this
	cdef Graph _G

	def __init__(self, *args, **kwargs):
		if type(self) == GraphLayoutAlgorithm:
			raise RuntimeError("Error, you may not use GraphLayoutAlgorithm directly, use a sub-class instead")

	def __dealloc__(self):
		self._G = None # just to be sure the graph is deleted

	def numEdgeCrossings(self):
		""" 
		numEdgeCrossings()
		
		Computes approximation (in parallel) of the Spanning Edge Centrality.
		
		Returns
		-------
		int
			Number of edge crossings.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.numEdgeCrossings()

	def getCoordinates(self):
		""" 
		getCoordinates()
		
		Computes approximation (in parallel) of the Spanning Edge Centrality.
		
		Returns
		-------
		list(tuple(float, float))
			List of coordinates for each node.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		cdef vector[vector[double]] pairCoord = vector[vector[double]]()
		cdef vector[double] prXd
		pointCoord = self._this.getCoordinates()
		num_dim = pointCoord[0].getDimensions() 
		if num_dim == 3:
			for pt in pointCoord:
				prXd = [pt[0], pt[1], pt[2]]
				pairCoord.push_back(prXd)
		elif num_dim == 2:
			for pt in pointCoord:
				prXd = [pt[0], pt[1]]
				pairCoord.push_back(prXd)
		else:
			raise RuntimeError("Currently graph layouting only supports 2D or 3D coordinates.")
		return pairCoord

	def run(self):
		"""
		run()

		Executes the graph layout algorithm.

		Returns
		-------
		networkit.viz.GraphLayoutAlgorithm
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		with nogil:
			self._this.run()
		return self

	def writeGraphToGML(self, path):
		"""
		writeGraphToGML(path)

		Writes the graph and its layout to a .gml file at the specified path.

		Parameters
		----------
		path: str
			Path where the graph file should be created.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.writeGraphToGML(stdstring(path))

	def writeKinemage(self, string path):
		"""
		writeKinemage(path)

		Writes the graph and its layout to a file at the specified path.

		Parameters
		----------
		path: str
			Path where the graph file should be created.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.writeKinemage(stdstring(path))

cdef extern from "<networkit/viz/MaxentStress.hpp>" namespace "NetworKit::MaxentStress":

	enum _GraphDistance "NetworKit::MaxentStress::GraphDistance":
		EDGE_WEIGHT,
		ALGEBRAIC_DISTANCE

class GraphDistance:
	EdgeWeight = EDGE_WEIGHT
	AlgebraicDistance = ALGEBRAIC_DISTANCE

cdef extern from "<networkit/viz/MaxentStress.hpp>" namespace "NetworKit::MaxentStress":

	enum _LinearSolverType "NetworKit::MaxentStress::LinearSolverType":
		LAMG,
		CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER,
		CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER

class LinearSolverType:
	Lamg = LAMG
	ConjugateGradientIdentityPreconditioner = CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER
	ConjugateGradientDiagonalPreconditioner = CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER

cdef extern from "<networkit/viz/MaxentStress.hpp>":

	cdef cppclass _MaxentStress "NetworKit::MaxentStress" (_GraphLayoutAlgorithm[double]):
		_MaxentStress(_Graph G, count dim, count k, double tolerance, _LinearSolverType linearSolverType, bool_t fastComputation, _GraphDistance graphDistance) except +
		_MaxentStress(_Graph G, count dim, const vector[Point[double]] coordinates, count k, double tolerance, _LinearSolverType linearSolverType, bool_t fastComputation, _GraphDistance graphDistance) except +
		void run() except +
		void scaleLayout() except +
		double computeScalingFactor() except +
		double fullStressMeasure() except +
		double maxentMeasure() except +
		double meanDistanceError() except +
		double ldme() except +
		void setQ(double q) except +
		void setAlpha(double alpha) except +
		void setAlphaReduction(double alphaReduction) except +
		void setFinalAlpha(double finalAlpha) except +
		void setConvergenceThreshold(double convThreshold) except +
		double getRhs() except +
		double getApproxEntropyTerm() except +
		double getSolveTime() except +


cdef class MaxentStress (GraphLayoutAlgorithm):
	"""
	MaxentStress(G, dim, k, coordinates=list(), tolerance=1e-5, linearSolverType=networkit.viz.LinearSolverType.Lamg, fastComputation=False, graphDistance=networkit.community.Normalization.EdgeWeight)

	Implementation of MaxentStress by Gansner et al. using a Laplacian system solver.
  	@see Gansner, Emden R., Yifan Hu, and Steve North. "A maxent-stress model for graph layout."
	Visualization and Computer Graphics, IEEE Transactions on 19, no. 6 (2013): 927-940.

	Parameter :code:`graphDistance` can be one of the following:

	- networkit.viz.GraphDistance.EdgeWeight
	- networkit.viz.GraphDistance.AlgebraicDistance

	Parameter :code:`linearSolverType` can be one of the following:

	- networkit.viz.LinearSolverType.Lamg
	- networkit.viz.LinearSolverType.ConjugateGradientIdentityPreconditioner
	- networkit.viz.LinearSolverType.ConjugateGradientDiagonalPreconditioner

	Parameters
	----------
	G : networkit.Graph
		The (connected) graph to be handled.
	dim: int
		Number of dimensions.
	k: int
		Node distance to take into account for computation. The higher k, the longer computation takes to complete.
	coordinates: list(tuple(float, float)), optional
		Fixed coordinates. Default: list()
	tolerance: float, optional
		The tolerance of the solver. Default: 1e-5
	linearSolverType: networkit.viz.LinearSolverType, optional
		The type of linear solver. Default: networkit.viz.LinearSolverType.Lamg
	fastComputation: bool, optional
		Decides whether or not slightly faster computation should be employed, leading to slightly worse results. Default: False
	graphDistance: networkit.viz.GraphDistance, optional
		Decides what type of graph distance should be utilised. Default: networkit.community.GraphDistance.EdgeWeight
	"""

	LAMG = 0
	CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER = 1
	CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER = 2
	EDGE_WEIGHT = 0
	ALGEBRAIC_DISTANCE = 1

	def __cinit__(self, Graph G, count dim, count k, vector[pair[double, double]] coordinates = [], double tolerance = 1e-5, _LinearSolverType linearSolverType = LAMG, bool_t fastComputation = False, _GraphDistance graphDistance = EDGE_WEIGHT):
		cdef Point[double] p = Point[double](0, 0)
		cdef vector[Point[double]] pointCoordinates = vector[Point[double]]()

		for pr in coordinates:
			p = Point[double](pr.first, pr.second)
			pointCoordinates.push_back(p)

		if (coordinates.size() != 0):
			self._this = new _MaxentStress(G._this, dim, pointCoordinates, k, tolerance, linearSolverType, fastComputation, graphDistance)
		else:
			self._this = new _MaxentStress(G._this, dim, k, tolerance, linearSolverType, fastComputation, graphDistance)

	def __dealloc__(self):
		del self._this

	def scaleLayout(self):
		"""
		scaleLayout()	

		Scale the layout computed by run() by a scalar s to minimize :math:`\sum_{u,v \in V} w_{uv} (s ||x_u - x_v|| - d_{uv}||)^2`.
		"""
		(<_MaxentStress*>(self._this)).scaleLayout()
		return self

	def computeScalingFactor(self):
		"""
		computeScalingFactor()

		Computes a scalar s s.t. :math:`\sum_{u,v \in V} w_{uv} (s ||x_u - x_v|| - d_{uv}||)^2` is minimized.
		"""
		return (<_MaxentStress*>(self._this)).computeScalingFactor()

	def fullStressMeasure(self):
		"""
		fullStressMeasure()

		Computes the full stress measure of the computed layout with run().
		"""
		return (<_MaxentStress*>(self._this)).fullStressMeasure()

	def maxentMeasure(self):
		"""
		maxentMeasure()
		
		Computes the maxent stress measure for the computed layout with run().
		"""
		return (<_MaxentStress*>(self._this)).maxentMeasure()

	def meanDistanceError(self):
		"""
		meanDistanceError()

		Computes mean distance error.
		"""
		return (<_MaxentStress*>(self._this)).meanDistanceError()

	def ldme(self):
		"""
		ldme()
		
		Computes the ldme.
		"""
		return (<_MaxentStress*>(self._this)).ldme()

	def setQ(self, double q):
		"""
		setQ(q)

		Set parameter q.

		Parameters
		----------
		q : float
			New parameter value.
		"""
		(<_MaxentStress*>(self._this)).setQ(q)
		return self

	def setAlpha(self, double alpha):
		"""
		setAlpha(alpha)

		Set parameter alpha.

		Parameters
		----------
		alpha : float
			New parameter value.
		"""
		(<_MaxentStress*>(self._this)).setAlpha(alpha)
		return self

	def setAlphaReduction(self, double alphaReduction):
		"""
		setAlphaReduction(alphaReduction)

		Set parameter alphaReduction.

		Parameters
		----------
		alphaReduction : float
			New parameter value.
		"""
		(<_MaxentStress*>(self._this)).setAlphaReduction(alphaReduction)
		return self

	def setFinalAlpha(self, double finalAlpha):
		"""
		setFinalAlpha(finalAlpha)

		Set parameter finalAlpha.

		Parameters
		----------
		finalAlpha : float
			New parameter value.
		"""
		(<_MaxentStress*>(self._this)).setFinalAlpha(finalAlpha)
		return self

	def setConvergenceThreshold(self, double convThreshold):
		"""
		setConvergenceThreshold(convThreshold)

		Set parameter convThreshold.

		Parameters
		----------
		convThreshold : float
			New parameter value.
		"""
		(<_MaxentStress*>(self._this)).setConvergenceThreshold(convThreshold)
		return self

	def getRhs(self):
		"""
		getRhs

		Returns rhs value.

		Returns
		-------
		float
			The parameter value.
		"""
		return (<_MaxentStress*>(self._this)).getRhs()

	def getApproxEntropyTerm(self):
		"""
		getApproxEntropyTerm()

		Returns entropy term value.

		Returns
		-------
		float
			The parameter value.
		"""
		return (<_MaxentStress*>(self._this)).getApproxEntropyTerm()

	def getSolveTime(self):
		"""
		getSolveTime

		Returns solve time value.

		Returns
		-------
		float
			The parameter value.
		"""
		return (<_MaxentStress*>(self._this)).getSolveTime()

cdef extern from "<networkit/viz/PivotMDS.hpp>":

	cdef cppclass _PivotMDS "NetworKit::PivotMDS" (_GraphLayoutAlgorithm[double]):
				_PivotMDS(_Graph G, count dim, count numberOfPivots) except +
				void run() except +


cdef class PivotMDS (GraphLayoutAlgorithm):
	"""
	PivotMDS(dim, numberOfPivots)

	Implementation of PivotMDS proposed by Brandes and Pich.

	Parameters
	----------
	G: networkit.Graph
		The graph to be handled by the algorithm.
	dim: int
		Number of dimensions.
	numberOfPivots: int
		Number of pivots for the algorithm.
	"""

	def __cinit__(self, Graph G, count dim, count numberOfPivots):
		self._this = new _PivotMDS(G._this, dim, numberOfPivots)

	def __dealloc__(self):
		del self._this
