import scipy.sparse

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

