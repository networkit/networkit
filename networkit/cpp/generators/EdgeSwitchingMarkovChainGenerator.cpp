// networkit-format

#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/randomization/EdgeSwitching.hpp>

namespace NetworKit {

EdgeSwitchingMarkovChainGenerator::EdgeSwitchingMarkovChainGenerator(
    const std::vector<count> &sequence, bool ignoreIfRealizable)
    : StaticDegreeSequenceGenerator(sequence), ignoreIfRealizable(ignoreIfRealizable) {}

Graph EdgeSwitchingMarkovChainGenerator::generate() {
    Graph result(HavelHakimiGenerator(seq, ignoreIfRealizable).generate());
    const count neededSwaps = result.numberOfEdges() * 10;

    EdgeSwitchingInPlace edgeSwitching(result, neededSwaps);
    edgeSwitching.run();

    return result;
}

} // namespace NetworKit
