""" This module deals with the conversion of graphs into matrices and linear algebra operations on graphs """


__author__ = "Christian Staudt"

# local imports

# external imports
try:
	import scipy.sparse
except ImportError:
	print("module 'scipy' not available -- some functionality will be restricted")
try:
	import numpy as np
except ImportError:
	print("module 'numpy' not available -- some functionality will be restricted")


def column(matrix, i):
	return [row[i] for row in matrix]


def adjacencyMatrix(G, matrixType="sparse"):
	""" Get the adjacency matrix of the graph `G`.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	matrixType : string
		represent"sparse" or "dense"

	Returns
	-------
	:py:class:`scipy.sparse.csr_matrix`
		The adjacency matrix of the graph.
	"""
	n = G.upperNodeIdBound()
	if matrixType is "sparse":
		A = scipy.sparse.lil_matrix((n,n))
	elif matrixType is "dense":
		A = np.zeros(shape=(n,n))
	else:
		raise InputError("unknown matrix type: '{0}'".format(matrixType))
	# TODO: replace .edges() with efficient iterations
	if G.isWeighted():
		if G.isDirected():
			def processEdge(u,v,w,id):
				A[u, v] = w
		else:
			def processEdge(u,v,w,id):
				A[u, v] = w
				A[v, u] = w
	else:
		if G.isDirected():
			def processEdge(u,v,w,id):
				A[u, v] = 1
		else:
			def processEdge(u,v,w,id):
				A[u, v] = 1
				A[v, u] = 1
	G.forEdges(processEdge)
	if matrixType is "sparse":
		A = A.tocsr()  # convert to CSR for more efficient arithmetic operations
	return A

def laplacianMatrix(G):
	""" Get the laplacian matrix of the graph `G`.

	Parameters
	----------
	G : networkit.Graph
		The graph.

	Returns
	-------
	lap : ndarray
		The N x N laplacian matrix of graph.
	diag : ndarray
		The length-N diagonal of the laplacian matrix.
		diag is returned only if return_diag is True.
	"""
	A = adjacencyMatrix(G)
	return scipy.sparse.csgraph.laplacian(A)

def PageRankMatrix(G, damp=0.85):
	"""
	Builds the PageRank matrix of the undirected Graph `G`. This matrix corresponds with the
	PageRank matrix used in the C++ backend.


	Parameters
	----------
	G : networkit.Graph
		The graph.
	damp:
		Damping factor of the PageRank algorithm (0.85 by default)
	Returns
	-------
	pr : ndarray
		 The N x N page rank matrix of graph.
	"""
	A = adjacencyMatrix(G)

	n = G.numberOfNodes()
	stochastify = scipy.sparse.lil_matrix((n,n))
	for v in G.nodes():
		neighbors = G.degree(v)
		if neighbors == 0:
			stochastify[v,v] = 0.0	# TODO: check correctness
		else:
			stochastify[v,v] = 1.0 / neighbors
	stochastify = stochastify.tocsr()

	stochastic = A * stochastify

	dampened = stochastic * damp

	teleport = scipy.sparse.identity(G.numberOfNodes(), format="csr") * ((1 - damp) / G.numberOfNodes())

	return dampened + teleport

def symmetricEigenvectors(matrix, cutoff=-1, reverse=False):
	"""
	Computes eigenvectors and -values of symmetric matrices.

	Parameters
	----------
	matrix : sparse matrix
			 The matrix to compute the eigenvectors of
	cutoff : int
			 The maximum (or minimum) magnitude of the eigenvectors needed
	reverse : boolean
			  If set to true, the smaller eigenvalues will be computed before the larger ones

	Returns
	-------
	pr : ( [ float ], [ ndarray ] )
		 A tuple of ordered lists, the first containing the eigenvalues in descending (ascending) magnitude, the
		 second one holding the corresponding eigenvectors

	"""
	if cutoff == -1:
		cutoff = matrix.shape[0] - 2

	if reverse:
		mode = "SA"
	else:
		mode = "LA"

	w, v = scipy.sparse.linalg.eigsh(matrix, cutoff + 1, which=mode)

	orderlist = zip(w, range(0, len(w)))
	orderlist = sorted(orderlist)

	orderedW = column(orderlist, 0)
	orderedV = [v[:,i] for i in column(orderlist, 1)]

	return (orderedW, orderedV)

def eigenvectors(matrix, cutoff=-1, reverse=False):
	"""
	Computes eigenvectors and -values of matrices.

	Parameters
	----------
	matrix : sparse matrix
			 The matrix to compute the eigenvectors of
	cutoff : int
			 The maximum (or minimum) magnitude of the eigenvectors needed
	reverse : boolean
			  If set to true, the smaller eigenvalues will be computed before the larger ones

	Returns
	-------
	pr : ( [ float ], [ ndarray ] )
		 A tuple of ordered lists, the first containing the eigenvalues in descending (ascending) magnitude, the
		 second one holding the corresponding eigenvectors

	"""
	if cutoff == -1:
		cutoff = matrix.shape[0] - 2

	if reverse:
		mode = "SR"
	else:
		mode = "LR"

	w, v = scipy.sparse.linalg.eigs(matrix, cutoff + 1, which=mode)

	orderlist = zip(w, range(0, len(w)))
	orderlist = sorted(orderlist)

	orderedW = column(orderlist, 0)
	orderedV = [v[:,i] for i in column(orderlist, 1)]

	return (orderedW, orderedV)

def laplacianEigenvectors(G, cutoff=-1, reverse=False):
	if G.isDirected():
		return eigenvectors(laplacianMatrix(G), cutoff=cutoff, reverse=reverse)
	else:
		return symmetricEigenvectors(laplacianMatrix(G), cutoff=cutoff, reverse=reverse)

def adjacencyEigenvectors(G, cutoff=-1, reverse=False):
	if G.isDirected():
		return eigenvectors(adjacencyMatrix(G), cutoff=cutoff, reverse=reverse)
	else:
		return symmetricEigenvectors(adjacencyMatrix(G), cutoff=cutoff, reverse=reverse)

def laplacianEigenvector(G, order, reverse=False):
	if G.isDirected():
		spectrum = eigenvectors(laplacianMatrix(G), cutoff=order, reverse=reverse)
	else:
		spectrum = symmetricEigenvectors(laplacianMatrix(G), cutoff=order, reverse=reverse)

	return (spectrum[0][order], spectrum[1][order])

def adjacencyEigenvector(G, order, reverse=False):
	if G.isDirected():
		spectrum = eigenvectors(adjacencyMatrix(G), cutoff=order, reverse=reverse)
	else:
		spectrum = symmetricEigenvectors(adjacencyMatrix(G), cutoff=order, reverse=reverse)

	return (spectrum[0][order], spectrum[1][order])
