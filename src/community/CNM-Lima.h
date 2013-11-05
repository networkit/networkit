/*
 * CNM.h
 *
 *  Created on: 25.02.2013
 *      Author: Fellipe Lima
 */

#ifndef CNM_H_
#define CNM_H_

#include "Clusterer.h"
//#include "../base/IndexMap.h"
//#include "../independentset/Luby.h"
//#include "../Globals.h"

//#include "omp.h"

//#include <algorithm> // max-heap
//#include <vector>
//#include <map> // fraction of edge ends for each single community
//#include <iostream>

namespace NetworKit {

	// forward declarations
	class ModularityDelta;

	/**
 	 * This class represents a row of a sparse matrix.
 	 */
	class Row
	{
		private:
			// index of row
			index m_index;
			// balanced binary tree of deltas belonging to the row
			// TODO: auf "set" umsetzen, so dass er richtige Baeume verwendet
			std::vector<ModularityDelta*> m_DeltaQs;

		public:
			Row(index idx);

			// accessors
			index GetIndex() const { return m_index; };
			void GetDelta(index col, double& delta);
			std::vector<ModularityDelta*> GetDeltaList() { return m_DeltaQs; };
			void GetMaxDelta(ModularityDelta& delta);
			ModularityDelta GetMaxDelta();

			// mutators
			void RemoveCommunityReferences(index idx);
			void AddDelta(ModularityDelta* delta);

			// aux
			void PrintColumns(index dim);
	};

	/**
	 * This class represents a modularity delta (DeltaQ) between two communities i and j.
	 */
	class ModularityDelta
	{
		private:
			// index of first community
			index m_i;
			// index of second community
			index m_j;
			// value of modularity delta (DeltaQ)
			double m_DeltaValue;

		public:
			ModularityDelta();
			ModularityDelta(index i, index j, count iDegree, count jDegree, count edgeCount);

			// mutators
			void SetRow(index i) { m_i = i; };
			void SetCol(index j) { m_j = j; };
			void SetDelta(double delta);

			// accessors
			double GetDeltaValue() const { return m_DeltaValue; };
			index GetRow() const { return m_i; };
			index GetCol() const { return m_j; };

			// used for building the heap
			struct CompareDelta
			{
				bool operator()(ModularityDelta* a, ModularityDelta* b) const
				{
					return a->GetDeltaValue() < b->GetDeltaValue();
				}
			};
	};

	class CNM: public NetworKit::Clusterer
	{
		private:
			// each row is a balanced binary tree, and the list of rows composes the sparse matrix
			std::vector<Row*> m_SparseMatrix;
			// max-heap for fast fetching of the maximum delta among all rows
			std::vector<ModularityDelta*> m_MaxHeap;
			// fraction of edge ends for each community
			std::vector<double> m_EdgeEndFraction;

		public:
			CNM();
			virtual ~CNM();
			virtual Clustering run(Graph& G);
			virtual void InitCNMData(Graph& G);

			// algorithm methods
			virtual void BuildClustering(Graph& G, Clustering& clustering);
			virtual void RemoveCommunityAndRebuildHeap(index removedCommunity, count n);
			virtual void UpdateData(Graph& G, Clustering clustering, node i, node j);
			virtual void UpdateMatrix(Graph& G, Clustering clustering, node i, node j);
			virtual void UpdateEdgeEndFractions(node i, node j);

			// aux functions
			bool IsConnected(Graph& G, Clustering clustering, node i, node k);
			virtual bool IsHeapValid();
			void PrintMatrix();
			virtual std::string toString() const;

			// accessors
			Row* GetRow(index idx);
			void FetchDelta(index row, index col, double& delta);
			void FetchFraction(index i, double& fraction);
	};

} /* namespace NetworKit */
#endif /* CNM_H_ */
