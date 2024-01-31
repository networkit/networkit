
#ifndef NETWORKIT_GRAPH_PARALLEL_GRAPH_BUILDER_HPP_
#define NETWORKIT_GRAPH_PARALLEL_GRAPH_BUILDER_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class ParallelGraphBuilder {
    struct HalfEdge {
        index source;
        index destination;
    };
    count n;
    bool directed;
    std::vector<std::vector<std::vector<HalfEdge>>> outEdgesPerThread;
    std::vector<std::vector<std::vector<HalfEdge>>> inEdgesPerThread;

public:
    ParallelGraphBuilder(count n, bool directed = false);

    void addEdgeParallel(index a, index b);

    void addHalfEdgeParallel(Unsafe, index source, index destination);

    void addHalfInEdgeParallel(Unsafe, index source, index destination);

    void addHalfOutEdgeParallel(Unsafe, index source, index destination);

    Graph toGraphParallel();

    void addEdgesToGraphParallel(Unsafe, Graph &G);
};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_PARALLEL_GRAPH_BUILDER_HPP_