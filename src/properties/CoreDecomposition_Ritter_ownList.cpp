/*
 * CoreDecomposition_Ritter_ownList.cpp
 */

#include "CoreDecomposition_Ritter_ownList.h"
#include <list>

namespace NetworKit {

std::vector<count> CoreDecomposition_Ritter_ownList::run(const Graph& G) {
	const count n = G.numberOfNodes();

	// save resulting coreness for every node
	std::vector<count> coreness(n, 0);

	// data structure A, array of double linked lists
	// the list a[i] should contain all nodes with degree i
	// maximum degree is n-1, minimum degree is 0
	ListEntry* a[n];
	for (index i = 0; i < n; i++)
		a[i] = NULL;
	
	// save degree for every node (is modified during the calculations)
	int degree[n];
	// every nodes gets one ListEntry, map from node to ListEntry to access the ListEntry of the node without searching the lists
	ListEntry* vptr[n];
	
	// iterate over all nodes in graph G and initialize a, degree and vptr
	G.forNodes([&](node v) {
		const int d = G.degree(v);

		// create list entry
		struct ListEntry* e = (struct ListEntry*) malloc(sizeof(struct ListEntry));
		e->v = v;
		prepend(a, d, e);
		
		degree[v] = d;
		vptr[v] = e;
	});

	index d_min = 0; // a[d_min] will be our first not empty list
	count i = 1; // current coreness

	// choose a node with minimal degree, remove it from a, save the coreness i and decrease the degree of all adjacent nodes
	// when d_min >= i, and a[i] empty (and because d_min >= i all a[<= i] are empty) we are done with coreness i and continue with i+1
	for (count nodes_left = n; nodes_left > 0; nodes_left--) {
		// update d_min and i so that a[d_min] is still the first not empty list
		while (a[d_min] == NULL && d_min < n) {
			d_min++;
		}
		if (d_min >= i) {
			// if we are done the coreness i we can jump to coreness d_min + 1,
			// because if d_min > i + 1 there cannot be nodes left with coreness < d_min + 1
			i = d_min + 1;
		}

		// the first node in a[d_min] is the node we are "removing" (and all it edges)
		const node v_del = a[d_min]->v;

		// delete node from list
		// free(pop(a, d_min));
		pop(a, d_min);
		vptr[v_del] = NULL;

		// mark node as done and save coreness
		degree[v_del] = 0;
		coreness[v_del] = i - 1;

		// decrease degree of all adjacent nodes
		// keep structure of a (all nodes of a[i] have degree i)
		G.forEdgesOf(v_del, [&] (node v, node w) { // v == v_del
			const index d_old = degree[w];
			if (d_old == 0) {
				// if the degree is 0 w was already removed from a
				return;
			}

			const index d_new = d_old - 1;
			degree[w] = d_new;

			// move w from a[d_old] to a[d_new]
			remove(a, d_old, vptr[w]);
			prepend(a, d_new, vptr[w]);

			if (d_new < d_min) {
				d_min = d_new;
			}
		});
	}

	return coreness;
}

/*
 * Private helper methods for own double linked list implementation.
 */

/*
 * Prints the data structure a with all the nodes for every list.
 */
void CoreDecomposition_Ritter_ownList::printLists(ListEntry** a, int n) {
	for (int i = 0; i < n; i++) {
		if (a[i] != NULL) {
			std::cout << "a[" << i << "]=[";

			ListEntry* e = a[i];
			while (e != NULL) {
				std::cout << e->v << ", ";
				e = e->next;
			}

			std::cout << "]\n";
		}
	}
}

/*
 * Prepends e to a[d].
 */
void CoreDecomposition_Ritter_ownList::prepend(ListEntry** a, int d, ListEntry* e) {
	e->prev = NULL;
	e->next = a[d]; // can be NULL if a[d] is empty
	if (a[d] != NULL) {
		a[d]->prev = e;
	}
	a[d] = e;
}

/*
 * Removes e from a, e must be in a[d].
 */
void CoreDecomposition_Ritter_ownList::remove(ListEntry** a, int d, ListEntry* e) {
	if (e->prev == NULL) {
		// e is first entry
		if (e->next != NULL) {
			// e is not last entry => next node becomes first entry
			e->next->prev = NULL;
		}
		a[d] = e->next; // NULL if e is last entry
	} else {
		// e is not first entry
		e->prev->next = e->next;
		if (e->next != NULL) {
			// e is not last entry
			e->next->prev = e->prev;
		}
	}
	e->prev = NULL;
	e->next = NULL;
}

/*
 * Removes and returns the first element from a[d]. a[d] must not be empty.
 */
CoreDecomposition_Ritter_ownList::ListEntry* CoreDecomposition_Ritter_ownList::pop(ListEntry** a, int d) {
	if (a[d] == NULL) {
		return NULL;
	}
	ListEntry* e = a[d];
	if (e->next != NULL) {
		// e isn't last entry, next entry becomes first element
		e->next->prev = NULL;
	}
	a[d] = e->next; // NULL if e is last element
	return e;
}

} /* namespace NetworKit */
