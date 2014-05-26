/*
 * ClusteringCuts.h
 *
 *  Created on: 23.10.2013
 *      Author: Weijian Ji (weijian.ji@student.kit.edu)
 */

#ifndef CLUSTERINGCUTS_H_
#define CLUSTERINGCUT_H_

#include "../graph/Graph.h"
#include "../graph/NodeMap.h"
#include "../structures/Partition.h"
#include <set>
#include <unordered_map>
#include <utility> // for including Pair

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}



namespace NetworKit {

typedef uint64_t count; // more expressive name for an integer quantity
typedef std::pair<node, node>  edge; // more expressive name for an edge
typedef std::set<std::pair<node, node> > cutEdges; // more expressive name for group of cutting edges
typedef std::unordered_map<std::pair<index, index>, cutEdges > cuttingMap; // map from cluster pair to their clustering cut edges
class ClusteringCut {

public:

	ClusteringCut();

	virtual ~ClusteringCut();

	/**
	 * create matrix map from cluster pair to their cutting edges
	 */
	cuttingMap getClusterToCutMatrixMap(const Graph& G, const Partition& zeta);

	/*
	 * get all cut edges of the clustering
	 */
	cutEdges getAllCutEdges(const Graph& G, const Partition& zeta);
	/*
	 * get cut between two clusters
	 */
	cutEdges getInterClusterCutEdges(index c1, index c2, const Graph& G, const Partition& zeta);

	/**
	 * return size of the clustering cutting edges between cluster of index1 and of index2
	 */
	count getInterClusterCutEdgesSize(index c1, index c2, const Graph& G, const Partition& zeta) const;

	/**
	 * return the total size of the clustering cutting edges
	 */
	count getTotalCutEdgesSize(const Graph& G, const Partition& zeta) const;

};
}/* namespace NetworKit */
#endif /* CLUSTERINGCUT_H_ */
