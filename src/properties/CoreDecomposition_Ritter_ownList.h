/*
 * CoreDecomposition_Ritter_ownList.h
 */

#ifndef COREDECOMPOSITION_RITTER_OWNLIST_H_
#define COREDECOMPOSITION_RITTER_OWNLIST_H_

#include <vector>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition_Ritter_ownList {
public:
	CoreDecomposition_Ritter_ownList();
	virtual ~CoreDecomposition_Ritter_ownList();

	/**
	 * @return k-core decomposition of graph @a G.
	 */
	std::vector<count> run(const Graph& G);

private:
	struct ListEntry
	{
		ListEntry* prev;
		ListEntry* next;
		node v;
	};

	void printLists(ListEntry** a, int n);
	void prepend(ListEntry** a, int d, ListEntry* e);
	void remove(ListEntry** a, int d, ListEntry* e);
	ListEntry* pop(ListEntry** a, int d);
};

} /* namespace NetworKit */
#endif /* COREDECOMPOSITION_RITTER_OWNLIST_H_ */
