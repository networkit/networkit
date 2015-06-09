/*
 * SpanningTree.h
 *
 *  Created on: 03 May 2014
 *      Author: dhoske
 */

#ifndef SPANNING_TREE_H
#define SPANNING_TREE_H

#include "RootedTree.h"

namespace NetworKit {
namespace SDD {

/** @defgroup sdd-st SpanningTree
 *  @ingroup sdd
 *  Algorithms for computing spanning trees (STs).
 *
 *  We provide implementations for a minimum distance STs
 * (@ref minDistanceST with Dijkstra), a minimum weight STs
 * (@ref minWeightST with Kruskal) and a low-stretch STs (@ref lowStretchST,
 * algorithm described in <a>&ldquo;Lower-Stretch Spanning Trees&rdquo;</a>
 * by Elkin et al.).
 *  @{ */

/**
 * Algorithm for computing a spanning tree of a connected graph.
 */
using SpanningTreeAlgo = RootedTree(*)(const Graph&, node);

/**
 * Returns the minimum distance spanning tree of the
 * connected graph @a G rooted at @a u in \f$O((m+n)\log(n))\f$O time
 * with Dijkstra's algorithm. @a G needs to be connected.
 */
RootedTree minDistanceST(const Graph& G, node u);

/**
 * Returns a minimum-weight spanning tree of @a G rooted
 * at @a u in \f$O(n\log(n))\f$ time with Kruskal's algorithm.
 * @a G needs to be connected.
 */
RootedTree minWeightST(const Graph& G, node u);

/**
 * Returns a spanning tree of @a G rooted at @a u that
 * has stretch in \f$O(m\log(n)\log\log(n))\f$ in
 * \f$O(m\log(n)^2\log\log(n))\f$ time. @a G needs to be connected.
 * This is an implementation an algorithm by
 * <a href="http://arxiv.org/abs/cs.DS/0411064">Elkin et al.</a>.
 *
 * @bug As noted in <a href="http://www.cs.elte.hu/blobs/diplomamunkak/bsc_alkmat/2014/papp_pal_andras.pdf">
 * Pál András Papp's</a> bachelor thesis, there is a small error in the proof of
 * the stretch bound that can be fixed. We did not implement the fixed version, since
 * it yields worse stretch in practice.
 */
RootedTree lowStretchST(const Graph& G, node u);

/**
 * Returns a special spanning tree on a 2D-grid that is easily shown to have
 * \f$O(n\log n)\f$ stretch. See the thesis for more details.
 *
 * @warning @a G is expected to a be a square grid and this is not checked.
 */
RootedTree specialGridST(const Graph& G, node u);

/** @} */

}
}

#endif
