import scipy.sparse
import numpy as np
import math
from _NetworKit import Partition, Graph

def adjacencyMatrix(G):
	""" Return the adjacency matrix of an undirected graph."""
	n = G.numberOfNodes()	
	
	'create adjacency matrix in coo format'	
	row = []
	col = []
	data = []

	for (u,v) in G.edges():		
		row.append(u)
		col.append(v)
		data.append(G.weight(u,v))

		# symmetric 		
		row.append(v)
		col.append(u)
		data.append(G.weight(v,u))

	return scipy.sparse.coo_matrix((data, (row,col)), shape=(n,n)).tocsr()

def laplacianMatrix(G):
	""" Return the laplacian matrix of an undirected graph"""
	A = adjacencyMatrix(G)
	return scipy.sparse.csgraph.laplacian(A)

def spectralPartitioning(G, k=2):
	partition = Partition(G.numberOfNodes()) # create new partition
	partition.setUpperBound(k)

	L = scipy.sparse.csc_matrix(laplacianMatrix(G)) # get the laplacian matrix of G in csc format
	#evals, evecs = scipy.sparse.linalg.eigsh(L, 2, which='SM') # Compute the two second smallest eigenvectors of L (eigenvalues must be >= 0, because L is always positive-semidefinite)
	#fiedlerVector = evecs[:,1] # take the second smallest eigenvector, the Fiedler vector
	#L = L.astype(np.float)
	X = scipy.rand(L.shape[0], 2)
	(evals,evecs) = scipy.sparse.linalg.lobpcg(L, X, M=None, tol=1e-12, largest=False, maxiter=500, verbosityLevel=0)
	fiedlerVector = evecs[:,1]


	median = sorted(fiedlerVector)[len(fiedlerVector) // 2 - 1] # compute median

	recursionDepth = math.log(k, 2)
	if (recursionDepth > 1):
		G1 = Graph(0, True)
		G2 = Graph(0, True)
		oldToNewIds1 = {}
		oldToNewIds2 = {}
		newToOldIds1 = {}
		newToOldIds2 = {}

		for idx in range(0, len(fiedlerVector)):			
			if (fiedlerVector[idx] <= median):
				addedNode = G1.addNode()
				oldToNewIds1[idx] = addedNode
				newToOldIds1[addedNode] = idx
			else:
				addedNode = G2.addNode()
				oldToNewIds2[idx] = addedNode
				newToOldIds2[addedNode] = idx

		for (u,v) in G.edges():
			if (u in oldToNewIds1):
				if (v in oldToNewIds1):
					newU = oldToNewIds1[u]
					newV = oldToNewIds1[v]
					G1.addEdge(newU, newV, G.weight(newU, newV))

			else:
				if (v in oldToNewIds2):
					newU = oldToNewIds2[u]
					newV = oldToNewIds2[v]
					G2.addEdge(newU, newV, G.weight(newU, newV))


		p1 = spectralPartitioning(G1, k // 2)
		p2 = spectralPartitioning(G2, k // 2)

		for i in range(0, k // 2):
			for partitionMember in p1.getMembers(i):
				partition.addToSubset(i, newToOldIds1[partitionMember])

			for partitionMember in p2.getMembers(i):
				partition.addToSubset(i + k // 2, newToOldIds2[partitionMember])

	else:
		for idx in range(0, len(fiedlerVector)): # fill the two partitions
		    if (fiedlerVector[idx] <= median):
		        partition.addToSubset(0, idx)
		    else:
		        partition.addToSubset(1, idx)

	
	return partition
