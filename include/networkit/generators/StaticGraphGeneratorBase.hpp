/*
 * StaticGraphGeneratorBase.hpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_STATIC_GRAPH_GENERATOR_BASE_HPP_
#define NETWORKIT_GENERATORS_STATIC_GRAPH_GENERATOR_BASE_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Abstract base class for static graph generators.
 */
template <typename GraphType>
class StaticGraphGeneratorBase {

public:
    /** Default destructor */
    virtual ~StaticGraphGeneratorBase() = default;

    virtual GraphType generate() = 0;
};

using StaticGraphGenerator = StaticGraphGeneratorBase<Graph>;
using StaticHypergraphGenerator = StaticGraphGeneratorBase<Hypergraph>;

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_STATIC_GRAPH_GENERATOR_BASE_HPP_
