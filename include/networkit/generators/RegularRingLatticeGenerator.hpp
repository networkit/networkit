/*
* RegularRingLatticeGenerator.h
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

#ifndef NETWORKIT_GENERATORS_REGULAR_RING_LATTICE_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_REGULAR_RING_LATTICE_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class RegularRingLatticeGenerator: public StaticGraphGenerator {

public:
    /**
    * Construct a undirected regular ring lattice.
    *
    * @param nNodes 		number of nodes in target graph
    * @param nNeighbors		number of neighbors on each side of a node
    */
    RegularRingLatticeGenerator(count nNodes, count nNeighbors);

    virtual Graph generate();

protected:
        count nNodes;
        count nNeighbors;

};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_REGULAR_RING_LATTICE_GENERATOR_HPP_
