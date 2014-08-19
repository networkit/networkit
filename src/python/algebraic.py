import scipy.sparse

def adjacencyMatrix(G):
	""" Return the adjacency matrix of an undirected graph."""
	n = G.numberOfNodes()
	A = scipy.sparse.lil_matrix((n,n))
	for (u, v) in G.edges():
		A[u, v] = G.weight(u, v)
		A[v, u] = G.weight(v, u)
	A = A.tocsr()
	return A

def laplacianMatrix(G):
	""" Return the laplacian matrix of an undirected graph"""
	A = adjacencyMatrix(G)
	return scipy.sparse.csgraph.laplacian(A)

