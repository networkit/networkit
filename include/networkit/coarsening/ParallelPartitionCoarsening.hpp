/*
 * ParallelPartitionCoarsening.h
 *
 *  Created on: 03.07.2014
 *      Author: cls
 */

#ifndef NETWORKIT_COARSENING_PARALLEL_PARTITION_COARSENING_HPP_
#define NETWORKIT_COARSENING_PARALLEL_PARTITION_COARSENING_HPP_

#include <networkit/Globals.hpp>
#include <networkit/coarsening/GraphCoarsening.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup coarsening
 */
class ParallelPartitionCoarsening: public GraphCoarsening {
public:
    ParallelPartitionCoarsening(const Graph& G, const Partition& zeta, bool useGraphBuilder = true);

    virtual void run();

private:
    const Partition& zeta;
    bool useGraphBuilder;
};

} /* namespace NetworKit */

#endif // NETWORKIT_COARSENING_PARALLEL_PARTITION_COARSENING_HPP_
