/*
 * CNM.cpp
 *
 *  Created on: Jun 10, 2013
 *      Author: fellipe
 */

#include "CNM-Lima.h"


namespace NetworKit {

#if 0

///<summary>
/// Constructor
///</summary>
///<param name = "index">Index of the row</param>
Row::Row(index idx) {
	m_index = idx;
}

///<summary>
/// Extract the modularity delta with maximum value within the row.
///</summary>
///<param name = "delta">Object where the maximum delta should be stored</param>
void Row::GetMaxDelta(ModularityDelta& delta) {
	double dMaxDeltaValue = 0.0;

	count iDeltaCount = m_DeltaQs.size();

	for (index i = 0; i < iDeltaCount; ++i) {
		ModularityDelta* currDelta = m_DeltaQs.at(i);

		double dCurrDeltaValue = currDelta->GetDeltaValue();

		if (dCurrDeltaValue > dMaxDeltaValue) {
			dMaxDeltaValue = dCurrDeltaValue;

			double dDeltaValue = currDelta->GetDeltaValue();
			index iRow = currDelta->GetRow();
			index iCol = currDelta->GetCol();

			delta.SetDelta(dDeltaValue);
			delta.SetRow(iRow);
			delta.SetCol(iCol);
		}
	}
}

///<summary>
/// Extract the modularity delta with maximum value within the row.
///</summary>
///<param name = "delta">Object where the maximum delta should be stored</param>
ModularityDelta Row::GetMaxDelta() {
	ModularityDelta* maxDelta = new ModularityDelta();
	double dMaxDeltaValue = 0.0;

	count iDeltaCount = m_DeltaQs.size();

	for (index i = 0; i < iDeltaCount; ++i) {
		ModularityDelta* currDelta = m_DeltaQs.at(i);

		double dCurrDeltaValue = currDelta->GetDeltaValue();

		if (dCurrDeltaValue > dMaxDeltaValue) {
			dMaxDeltaValue = dCurrDeltaValue;

			double dDeltaValue = currDelta->GetDeltaValue();
			index iRow = currDelta->GetRow();
			index iCol = currDelta->GetCol();

			maxDelta->SetDelta(dDeltaValue);
			maxDelta->SetRow(iRow);
			maxDelta->SetCol(iCol);
		}
	}

	return *maxDelta;
}

///<summary>
/// Get modularity delta for a position within the row.
///</summary>
///<param name = "col">Column index within the row, identifying the wished position</param>
///<param name = "delta">Object where the  delta should be stored</param>
void Row::GetDelta(index col, double& delta) {
	count iColumnCount = m_DeltaQs.size();

	delta = 0.0;

	for (index i = 0; i < iColumnCount; ++i) {
		ModularityDelta* currDelta = m_DeltaQs.at(i);
		index iCurrColumnI = currDelta->GetRow();
		index iCurrColumnJ = currDelta->GetCol();

		if (iCurrColumnI == col || iCurrColumnJ == col) {
			delta = currDelta->GetDeltaValue();
		}
	}
}

///<summary>
/// Add a new delta to the row.
///</summary>
///<param name = "delta">Delta to be added to the row</param>
void Row::AddDelta(ModularityDelta* delta) {
	if (delta != NULL) {
		m_DeltaQs.push_back(delta);
	}
}

///<summary>
/// Remove all references of a community.
///</summary>
///<param name = "index>Index of the community to have the references removed</param>
void Row::RemoveCommunityReferences(index idx) {
	count iColumnCount = m_DeltaQs.size();
	std::vector<ModularityDelta*> newRow;

	for (index i = 0; i < iColumnCount; ++i) {
		ModularityDelta* delta = m_DeltaQs.at(i);
		index iCurrI = delta->GetRow();
		index iCurrJ = delta->GetCol();
		double dCurrDelta = delta->GetDeltaValue();

		if (iCurrI != idx && iCurrJ != idx) {
			// generate a new modularity delta with the old data
			ModularityDelta* newDelta = new ModularityDelta();

			newDelta->SetDelta(dCurrDelta);
			newDelta->SetRow(iCurrI);
			newDelta->SetCol(iCurrJ);

			newRow.push_back(newDelta);
		}
	}

	m_DeltaQs.clear();

	for (index i = 0; i < newRow.size(); ++i) {
		m_DeltaQs.push_back(newRow.at(i));
	}
}

///<summary>
/// Constructor
///</summary>
ModularityDelta::ModularityDelta() {
}

///<summary>
/// Constructor for initialization step, where the clusters are vertices.
///</summary>
///<param name = "i">First vertex/community (row index of the sparse matrix)</param>
///<param name = "j">Second vertex/community (row index of the sparse matrix)</param>
///<param name = "ki">Unweighted degree of the first vertex/community</param>
///<param name = "kj">Unweighted degree of the second vertex/community</param>
///<param name = "m">Total number of edges of the graph</param>
ModularityDelta::ModularityDelta(index i, index j, count ki, count kj, count m) {
	m_i = i;
	m_j = j;
	m_DeltaValue = 0.0;

	if (ki >= 0 && kj >= 0) {
		// e(i,j)
		double dIntraCluster = (double) (1.0 / (2.0 * m));

		// p(i,j) = k(i)*k(j) / ((2*m)^2)
		double dExpectedValue = (double) ((ki * kj) / ((2.0 * m) * (2.0 * m)));

		// DeltaQ(i,j)
		m_DeltaValue = dIntraCluster - dExpectedValue;
	}
}

///<summary>
/// Set the delta value.
///</summary>
///<param name = "delta">Modularity delta</param>
void ModularityDelta::SetDelta(double delta) {
	m_DeltaValue = delta;
}

///<summary>
/// Constructor
///</summary>
CNM::CNM() {
}

///<summary>
/// Destructor
///</summary>
CNM::~CNM() {
}

///<summary>
/// Execute the CNM method for a graph.
///</summary>
///<param name = "G">Graph</param>
Clustering CNM::run(Graph& G) {
	// build initial clustering: each vertex lies within its own cluster
	count iGraphSize = G.numberOfNodes();
	Clustering result(iGraphSize);
	result.allToSingletons();

	// initialize data for CNM and build the agglomerative clustering
	InitCNMData(G);
	BuildClustering(G, result);

	return result;
}

///<summary>
/// This function calculates the (1) initial delta modularity values between every two communities - deltaQ(i,j).
/// In this case every vertex has its own cluster, so the problem is reduced to calculate the modularity between vertices.
/// Afterwards it calculates the (2) fraction of edge ends within each community - a(i). Furthermore it picks the (3) largest
/// deltaQ's of each single row and inserts it into a max-heap.
///</summary>
///<param name = "G">Graph</param>
void CNM::InitCNMData(Graph& G) {
	count n = G.numberOfNodes();
	count m = G.numberOfEdges(); // XXX: fixed by HM from numberOfNodes

	// 1) calculate the DeltaQ's between community i the neighborhood of communities j
	// and gather the maximum for each row
	for (index i = 0; i < n; ++i) {
		double dMaxDeltaValue = 0.0;
		ModularityDelta* maxDeltaInRow = NULL;
		Row* newRow = new Row(i);

		// calculate DeltaQ's
		for (index j = 0; j < n; ++j) {
			if (G.hasEdge(i, j)) {
				// degrees of the vertices for initial modularity delta calculation
				count ki = G.degree(i);
				count kj = G.degree(j);

				ModularityDelta* delta = new ModularityDelta(i, j, ki, kj, m);
				newRow->AddDelta(delta);

				// actualize the maximum delta lying in the row
				double dCurrDeltaValue = delta->GetDeltaValue();

				if (dCurrDeltaValue > dMaxDeltaValue) {
					dMaxDeltaValue = dCurrDeltaValue;
					maxDeltaInRow = delta;
				}
			}
		}

		// 2) set the initial fraction of edge ends attached to community i
		double dFraction = (G.degree(i) / (2.0 * m));
		m_EdgeEndFraction.push_back(dFraction);

		// 3) populate the max-heap and add the row to the sparse matrix:
		// the actual heap build occurs not before the complete traversing of all edges
		m_MaxHeap.push_back(maxDeltaInRow);
		m_SparseMatrix.push_back(newRow);
	}

	// 4) build the heap
	std::make_heap(m_MaxHeap.begin(), m_MaxHeap.end(),
			ModularityDelta::CompareDelta());
}

///<summary>
/// Build the agglomerative clustering with the data generated from the graph (G), firstly known from the function run.
///</summary>
///<param name = "G">Graph, for which the clustering will be built</param>
///<param name = "clustering">Clustering</param>
void CNM::BuildClustering(Graph& G, Clustering& clustering) {
	count iNodeCount = clustering.numberOfClusters();
	count iIterationCount = 0;

	count joinedCommunities[iNodeCount];

	do {
		if (IsHeapValid()) {
			// extract the community pair with largest modularity delta
			ModularityDelta* maxDelta = m_MaxHeap.front();
			if (maxDelta != NULL) {
				std::pop_heap(m_MaxHeap.begin(), m_MaxHeap.end());
				m_MaxHeap.pop_back();

				// join a community pair
				index i = maxDelta->GetRow();
				index j = maxDelta->GetCol();

				clustering.mergeClusters(i, j); // XXX changed from mergeClusters2 by HM
				joinedCommunities[iIterationCount] = i;

				// update the CNM data after joining communities i and j
				// the first community disappears and all nodes of it goes down into the second one
				UpdateData(G, clustering, i, j);
			}

			iIterationCount++;
		}
	} while (iIterationCount < iNodeCount);
}

///<summary>
/// Determines whether the heap property is satisfied for the max-heap.
///</summary>
bool CNM::IsHeapValid() {
	return std::is_heap(m_MaxHeap.begin(), m_MaxHeap.end(),
			ModularityDelta::CompareDelta());
}

///<summary>
/// Determine wether two communities are connected, i.e. whether two nodes belonging to both communities,
/// but different ones, are connected.
///</summary>
///<param name = "G">Graph</param>
///<param name = "Clustering">Clustering, used to gather all nodes of a community</param>
///<param name = "communityA">First community</param>
///<param name = "communityB">Second community</param>
bool CNM::IsConnected(Graph& G, Clustering clustering, index communityA,
		index communityB) {
	if (communityA == communityB) {
		return false;
	}

	std::vector<index> communityANodes = clustering.GetNodesForCluster(
			communityA);
	std::vector<index> communityBNodes = clustering.GetNodesForCluster(
			communityB);

	for (index i = 0; i < communityANodes.size(); ++i) {
		index noteFromA = communityANodes.at(i);

		for (index j = 0; j < communityBNodes.size(); ++j) {
			index noteFromB = communityBNodes.at(j);

			if (G.hasEdge(noteFromA, noteFromB)) {
				return true;
			}
		}

	}

	return false;
}

///<summary>
/// Get business logic object of row indexed index.
///</summary>
///<param name = "index">Index of the row to be fetched</param>
Row* CNM::GetRow(index idx) {
	Row* result = NULL;
	count iRowCount = m_SparseMatrix.size();
	for (index i = 0; i < iRowCount; ++i) {
		Row* currRow = m_SparseMatrix.at(i);
		index iCurrRowIndex = currRow->GetIndex();
		if (idx == iCurrRowIndex) {
			result = currRow;
		}
	}

	return result;
}

///<summary>
/// Remove community from the business logic and rebuild the max-heap containing the maxima of each row.
/// Firstly the row corresponding to the community will be removed, afterwards all references to this community occurring
/// in further rows are also removed.
///</summary>
///<param name = "removedCommunity">Index of the removed community, for which the CNM data will be cleaned</param>
///<param name = "nodeColunt">Total number of nodes</param>
void CNM::RemoveCommunityAndRebuildHeap(index removedCommunity, count nodeCount) {
	// remove the row
	// TODO: substitute this gadget for the actual "erase" of STL, which causes run time errors (???)
	count iRowCount = m_SparseMatrix.size();

	std::vector<Row*> newRows;

	for (index row = 0; row < iRowCount; ++row) {
		Row* currRow = m_SparseMatrix.at(row);

		if (currRow->GetIndex() != removedCommunity) {
			// add filtered row to the new (temporary) matrix
			newRows.push_back(currRow);

			// clean the columns of the row
			currRow->RemoveCommunityReferences(removedCommunity);
		}
	}

	// fill sparse matrix with filtered data
	m_SparseMatrix.clear();

	for (index i = 0; i < newRows.size(); ++i) {
		m_SparseMatrix.push_back(newRows.at(i));
	}

	// rebuild the heap with new row maximum data
	m_MaxHeap.clear();

	for (index i = 0; i < m_SparseMatrix.size(); ++i) {
		Row* currRow = m_SparseMatrix.at(i);

		// insert a row-maximum to the max-heap
		//ModularityDelta* maxRowDelta = new ModularityDelta();
		//currRow->GetMaxDelta(*maxRowDelta);
		ModularityDelta* maxRowDelta = new ModularityDelta();
		*maxRowDelta = currRow->GetMaxDelta();

		if (maxRowDelta != NULL) {
			index deltaCol = maxRowDelta->GetCol();
			index deltaRow = maxRowDelta->GetRow();

			if (deltaCol >= 0 && deltaRow >= 0 && deltaCol <= nodeCount
					&& deltaRow <= nodeCount) {
				//if (deltaCol != deltaRow)
				//{
				m_MaxHeap.push_back(maxRowDelta);
				//}
			}
		}
	}

	std::make_heap(m_MaxHeap.begin(), m_MaxHeap.end(),
			ModularityDelta::CompareDelta());
}

///<summary>
/// Print all the columns of a row.
///</summary>
///<param name = "dim">Dimension of the square sparse matrix</param>
void Row::PrintColumns(index dim) {
	std::vector<double> columnElements;

	for (index i = 0; i < dim; ++i) {
		columnElements.push_back(0.0);
	}

	if (dim >= m_DeltaQs.size()) {
		for (index i = 0; i < m_DeltaQs.size(); ++i) {
			index colIndex = m_DeltaQs.at(i)->GetCol();

			columnElements.insert(columnElements.begin() + colIndex,
					m_DeltaQs.at(i)->GetDeltaValue());
		}
	}

	for (index i = 0; i < dim; ++i) {
		std::cout << columnElements.at(i) << " ";
	}
}

///<summary>
/// Print the whole sparse matrix.
///</summary>
void CNM::PrintMatrix() {
	count rowCount = m_SparseMatrix.size();

	for (index row = 0; row < rowCount; ++row) {
		Row* currRow = m_SparseMatrix.at(row);
		currRow->PrintColumns(rowCount);
		std::cout << std::endl;
	}
}

///<summary>
/// Update the modularity delta.
///</summary>
///<param name = "G">Graph</param>
///<param name = "clustering">Clustering</param>
///<param name = "i">Index of the joining community</param>
///<param name = "j">Index of the joined community</param>
void CNM::UpdateMatrix(Graph& G, Clustering clustering, index i, index j) {
	//Row of the community, to which the community i joined
	Row* row = GetRow(j);

	if (row != NULL) {
		std::vector<ModularityDelta*> rowDeltaList = row->GetDeltaList();

		count iColumnCount = rowDeltaList.size();

		for (index columnIdx = 0; columnIdx < iColumnCount; ++columnIdx) {
			ModularityDelta* currDeltaQjk = rowDeltaList.at(columnIdx);

			index k = currDeltaQjk->GetCol();

			bool bConnection_ki = CNM::IsConnected(G, clustering, k, i);
			bool bConnection_kj = CNM::IsConnected(G, clustering, k, j);

			if (bConnection_ki && bConnection_kj) {
				double dNewDeltaQjk;
				double dDeltaQik;
				double dDeltaQjk;

				FetchDelta(i, k, dDeltaQik);
				FetchDelta(j, k, dDeltaQjk);

				dNewDeltaQjk = dDeltaQik + dDeltaQjk;
				currDeltaQjk->SetDelta(dNewDeltaQjk);
			} else if (bConnection_ki && !bConnection_kj) {
				double dNewDeltaQjk;
				double dDeltaQik;
				double dFraction_j;
				double dFraction_k;

				FetchDelta(i, k, dDeltaQik);
				FetchFraction(j, dFraction_j);
				FetchFraction(k, dFraction_k);

				dNewDeltaQjk = dDeltaQik - (2.0 * dFraction_j * dFraction_k);
				currDeltaQjk->SetDelta(dNewDeltaQjk);
			} else if (!bConnection_ki && bConnection_kj) {
				double dNewDeltaQjk;
				double dDeltaQjk;
				double dFraction_i;
				double dFraction_k;

				FetchDelta(j, k, dDeltaQjk);
				FetchFraction(i, dFraction_i);
				FetchFraction(k, dFraction_k);

				dNewDeltaQjk = dDeltaQjk - (2.0 * dFraction_i * dFraction_k);
				currDeltaQjk->SetDelta(dNewDeltaQjk);
			}
		}
	}
}

///<summary>
/// This function updates the whole CNM information after that communities i and j are joined.
/// In this case community i disappears and is joined to community j.
///</summary>
///<param name = "G">Graph</param>
///<param name = "clustering">Clustering</param>
///<param name = "i">Index of the joining community</param>
///<param name = "j">Index of the joined community</param>
void CNM::UpdateData(Graph& G, Clustering clustering, index i, index j) {
	count n = G.numberOfNodes();

	UpdateMatrix(G, clustering, i, j);
	RemoveCommunityAndRebuildHeap(i, n);
	UpdateEdgeEndFractions(i, j);
}

///<summary>
/// Update the fraction of edge end belonging to changed communities.
///</summary>
///<param name = "i">Index of the joining community</param>
///<param name = "j">Index of the joined community</param>
void CNM::UpdateEdgeEndFractions(index i, index j) {
	m_EdgeEndFraction[j] = m_EdgeEndFraction[j] + m_EdgeEndFraction[i];
	m_EdgeEndFraction[i] = 0.0;
}

///<summary>
/// Get the delta from the sparse matrix.
///</summary>
///<param name = "row">Index row</param>
///<param name = "col">Index column</param>
///<param name = "delta">Delta</param>
void CNM::FetchDelta(index row, index col, double& delta) {
	delta = 0.0;

	if (row != col) {
		count iRowCount = m_SparseMatrix.size();

		for (index i = 0; i < iRowCount; ++i) {
			Row* currRow = m_SparseMatrix.at(i);
			index currRowIndex = currRow->GetIndex();

			if (currRowIndex == row) {
				currRow->GetDelta(col, delta);
			}
		}

	}
}

///<summary>
/// Get the edge end fraction for a community.
///</summary>
///<param name = "i">Community</param>
///<param name = "fraction">Edge end fraction</param>
void CNM::FetchFraction(index i, double& fraction) {
	count iRowCount = m_EdgeEndFraction.size();

	for (index a = 0; a < iRowCount; ++a) {
		if (a == i) {
			fraction = m_EdgeEndFraction.at(i);
			return;
		}
	}
}

///<summary>
/// ?
///</summary>
std::string CNM::toString() const {
	return "CNM modularity clustering";
}

#endif

} /* namespace NetworKit */
