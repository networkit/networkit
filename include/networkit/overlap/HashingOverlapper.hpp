/*
 * HashingOverlapper.h
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_OVERLAP_HASHING_OVERLAPPER_HPP_
#define NETWORKIT_OVERLAP_HASHING_OVERLAPPER_HPP_

#include <functional>

#include <networkit/overlap/Overlapper.hpp>

namespace NetworKit {

/**
 * @ingroup overlap
 * Determines the overlap of multiple partitions by hashing partition identifiers.
 */
class HashingOverlapper: public Overlapper {

public:

    virtual Partition run(const Graph& G, const std::vector<Partition>& clusterings);

};

} /* namespace NetworKit */
#endif // NETWORKIT_OVERLAP_HASHING_OVERLAPPER_HPP_
