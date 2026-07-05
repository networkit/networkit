/*  GraphConceptsGTest.hpp
 *
 *  Created on: 05.07.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <string>
#include <gtest/gtest.h>
#include <networkit/Globals.hpp>
#include <networkit/graph/GraphConcepts.hpp>

namespace NetworKit {

TEST(GraphConceptsGTest, testGraphNodeConstraints) {
    static_assert(GraphNode<node>);
    static_assert(GraphNode<int>);
    static_assert(GraphNode<uint32_t>);
    static_assert(GraphNode<bool>);
    static_assert(!GraphNode<double>);

    SUCCEED();
}

TEST(GraphConceptsGTest, testGraphEdgeWeightConstraints) {
    static_assert(GraphEdgeWeight<edgeweight>);
    static_assert(GraphEdgeWeight<int>);
    static_assert(GraphEdgeWeight<double>);
    static_assert(GraphEdgeWeight<bool>);
    static_assert(!GraphEdgeWeight<std::string>);

    SUCCEED();
}

} // namespace NetworKit
