import scipy.sparse

from toolbox import column

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

def symmetricEigenvectors(matrix, cutoff=-1, reverse=False):
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
