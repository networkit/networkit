#include <networkit/subgraphs/GraphletsCounter.hpp>

#include "GraphletsCounterImpl.hpp"

namespace NetworKit {

using namespace NetworKit::SubgraphsDetails;

GraphletsCounter::GraphletsCounter(const Graph& G, unsigned K):
        G{&G}, k{K} {
}

void GraphletsCounter::run() {
    if(G->isDirected())
        throw std::runtime_error(
            "Graphlets counting is only implemented for undirected graphs"
        );
    // dispatch on the appropriate specialisation of GraphletsCounterImpl
    switch(k) {
        case 2u:
            counts = GraphletsCounterImpl<GraphletsSize::TWO>(G).getCounts();
            break;
        case 3u:
            counts = GraphletsCounterImpl<GraphletsSize::THREE>(G).getCounts();
            break;
        case 4u:
            counts = GraphletsCounterImpl<GraphletsSize::FOUR>(G).getCounts();
            break;
        default:
            throw std::runtime_error(
                "No algorithm is provided for " + std::to_string(k) +
                "-graphlets counting"
            );
    }
}
}  // namespace NetworKit
