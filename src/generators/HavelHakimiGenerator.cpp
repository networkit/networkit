/*
 * HavelHakimiGenerator.cpp
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include "HavelHakimiGenerator.h"

namespace NetworKit {

HavelHakimiGenerator::HavelHakimiGenerator(const std::vector<unsigned int>& sequence, bool skipTest) :
		seq(sequence) {
	std::sort(seq.begin(), seq.end(), std::greater<unsigned int>());
	if (skipTest) {
		realizable = true;
	}
	else {
		realizable = isRealizable();
	}
}

bool HavelHakimiGenerator::isRealizable() {
	DEBUG("check if sequence is realizable");
	realizable = true;
	count n = seq.size();

	/* First inequality. */
	unsigned int deg_sum = 0;
	for (count i = 0; i < n; ++i) {
		if (seq[i] < 0 || seq[i] >= n) {
			realizable = false;
			DEBUG("not realizable: ", seq[i], ", n: ", n);
			break;
		}
		deg_sum += seq[i];
	}

	if (deg_sum % 2 != 0) {
		DEBUG("not realizable");
		realizable = false;
	}

	/* Second inequality. */
	deg_sum = 0;
	for (unsigned int j = 0; j < n; ++j) {
		deg_sum += seq[j];

		/* sum of min(deg(i), j) for i from j + 1 to n - 1. */
		unsigned int min_deg_sum = 0;
		for (unsigned int i = j + 1; i < n; ++i) {
			min_deg_sum += std::min(seq[i], j + 1);
		}

		if (deg_sum > (j + 1) * j + min_deg_sum) {
			DEBUG("not realizable");
			realizable = false;
			break;
		}
	}

	return realizable;
}

bool HavelHakimiGenerator::getRealizable() const {
	return realizable;
}

Graph HavelHakimiGenerator::generate() {
	count n = seq.size();

	if (!realizable) {
		ERROR("Degree sequence not realizable or not checked for realizability yet");
		throw std::runtime_error("degree sequence not realizable or not checked for realizability");
	} else {
		DEBUG("Degree sequence is realizable, continue with generation algorithm.");
		Graph G(n);
		count numDegVals = (* std::max_element(seq.begin(), seq.end())) + 1;

		// Havel-Hakimi algorithm, adapted with appropriate data structure for linear time
		// bucket data structure: vector of lists
		// a list contains all vertices with the deficit corresponding to that list
		// in each iteration the first element of the highest non-empty list is chosen and removed
		// this vertex is connected to the following elements
		// when done, the other connected vertices are moved to the next list in reverse order

		typedef std::pair<count, node> DeficitAndNode;
		typedef std::vector<std::list<DeficitAndNode> > Buckets;

		// put nodes in appropriate lists
		Buckets nodesByDeficit(n);
		for(node v = 0; v < n; v++) {
			nodesByDeficit[seq[v]].push_front(std::make_pair(seq[v], v));
		}

		index maxDeficit = numDegVals - 1;
		while (maxDeficit > 0) {
//			DEBUG("maxDeficit: ", maxDeficit);

			// process node in largest bucket
			while(! nodesByDeficit[maxDeficit].empty()) {
				// get element
				std::list<DeficitAndNode>::iterator listIter = nodesByDeficit[maxDeficit].begin();
				count deficit = listIter->first;
				node currentVertex = listIter->second;

				// delete it
				nodesByDeficit[maxDeficit].pop_front();

				// connect corresponding vertex with the following ones
				index currentNeighborList = maxDeficit;
				std::stack<count> numToMove;

				while (deficit > 0) {
//					DEBUG("deficit: ", deficit);
//					DEBUG("currentNeighborList: ", currentNeighborList);
					count numDeleteFromCurrentList = 0;

					// search for candidates in current list
					for (auto elem : nodesByDeficit[currentNeighborList]) {
						// connect
						node nextNeighbor = elem.second;
//						DEBUG("add edge ", currentVertex, "-", nextNeighbor);
						G.addEdge(currentVertex, nextNeighbor);

						--deficit;
						++numDeleteFromCurrentList;

						if (deficit == 0) {
							++currentNeighborList; // due to -- a few lines below
							break;
						}
					}
					numToMove.push(numDeleteFromCurrentList);
					--currentNeighborList;
				}

				while (! numToMove.empty()) {
					// get head element
					count num = numToMove.top();
					numToMove.pop();

					// move this many items from current list to the next one
					for (index i = 0; i < num; ++i) {
						DeficitAndNode dan = nodesByDeficit[currentNeighborList].front();
						dan.first--;
						nodesByDeficit[currentNeighborList - 1].push_front(dan);
						nodesByDeficit[currentNeighborList].pop_front();
					}

					++currentNeighborList;
				}
//				DEBUG("adapt nodes in set");
			}
			maxDeficit--;
		}

		return G;
	}

}

} /* namespace NetworKit */
