# local imports

# external imports
import scipy.sparse
import numpy as np




def column(matrix, i):
    return [row[i] for row in matrix]


def adjacencyMatrix(G):
	""" Get the adjacency matrix of the undirected graph `G`.

	Parameters
	----------
	G : Graph
		The graph.

	Returns
	-------
	:py:class:`scipy.sparse.csr_matrix`
		The adjacency matrix of the graph.
	"""
	n = G.numberOfNodes()
	A = scipy.sparse.lil_matrix((n,n))
	for (u, v) in G.edges():
		A[u, v] = G.weight(u, v)
		A[v, u] = G.weight(v, u)
	A = A.tocsr()
	return A

def laplacianMatrix(G):
	""" Get the laplacian matrix of the undirected graph `G`.

	Parameters
	----------
	G : Graph
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

def PageRankMatrix(G, damp):
	"""
	Builds the PageRank matrix of the undirected Graph `G`. This matrix corresponds with the
	PageRank matrix used in the C++ backend.


	Parameters
	----------
	G : Graph
		The graph.

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

def laplacianEigenvectors(G, cutoff=-1, reverse=False):
	return symmetricEigenvectors(laplacianMatrix(G), cutoff=cutoff, reverse=reverse)

def adjacencyEigenvectors(G, cutoff=-1, reverse=False):
	return symmetricEigenvectors(adjacencyMatrix(G), cutoff=cutoff, reverse=reverse)

def laplacianEigenvector(G, order, reverse=False):
	spectrum = symmetricEigenvectors(laplacianMatrix(G), cutoff=order, reverse=reverse)

	return (spectrum[0][order], spectrum[1][order])

def adjacencyEigenvector(G, order, reverse=False):
	spectrum = symmetricEigenvectors(adjacencyMatrix(G), cutoff=order, reverse=reverse)

	return (spectrum[0][order], spectrum[1][order])
