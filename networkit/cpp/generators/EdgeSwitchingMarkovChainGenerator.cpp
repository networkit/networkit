#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>

NetworKit::EdgeSwitchingMarkovChainGenerator::EdgeSwitchingMarkovChainGenerator(const std::vector< NetworKit::count > &sequence, bool ignoreIfRealizable): StaticDegreeSequenceGenerator(sequence), ignoreIfRealizable(ignoreIfRealizable) {}

NetworKit::Graph NetworKit::EdgeSwitchingMarkovChainGenerator::generate() {
    Graph result(HavelHakimiGenerator(seq, ignoreIfRealizable).generate());

    count neededSwaps = result.numberOfEdges() * 10;
    count maxTry = neededSwaps * 2;
    count performedSwaps = 0;

    std::vector<node> nodeSelection;
    nodeSelection.reserve(result.numberOfEdges() * 2);

    result.forNodes([&](node u) {
        for (count i = 0; i < result.degree(u); ++i) {
            nodeSelection.push_back(u);
        }
    });

    for (count attempts = 0; attempts < maxTry && performedSwaps < neededSwaps; ++attempts) {
        node s1 = Aux::Random::choice(nodeSelection);
        node s2 = Aux::Random::choice(nodeSelection);

        if (s1 == s2) continue;

        node t1 = GraphTools::randomNeighbor(result, s1);
        node t2 = GraphTools::randomNeighbor(result, s2);

        if (t1 == t2 || s1 == t2 || s2 == t1) continue;

        if (result.hasEdge(s1, t2) || result.hasEdge(s2, t1)) continue; // FIXME make efficient!

        result.swapEdge(s1, t1, s2, t2);

        ++performedSwaps;
    }

    if (performedSwaps < neededSwaps) {
        INFO("Did only perform ", performedSwaps, " instead of ", neededSwaps, " edge swaps but made ", maxTry, " attempts to swap an edge");
    }

    return result;
}
