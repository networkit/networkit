#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ConfigurationModel.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/graph/Graph.hpp>

#include <chrono>
#include <iostream>
#include <thread>

namespace NetworKit {

ConfigurationModel::ConfigurationModel(const std::vector<count> &sequence)
    : StaticDegreeSequenceGenerator(sequence) {

    std::sort(seq.begin(), seq.end(), std::greater<count>());

    // Throw an error if the sequence isn't realizable:
    if (!isRealizable()) {
        throw std::runtime_error("Degree sequence is not realizable");
    }
}

Graph ConfigurationModel::generate() {
    do {
        Graph result(seq.size());
        std::vector<index> urn{};
        std::vector<count> remaining(seq.size());
        std::vector<bool> alreadyTaken(seq.size(), false);
        index remainingNodes = seq.size();

        // Fill the urn, s.t., the node of any degree is represented exactly degree times in the urn
        for (index i = 0; i < seq.size(); ++i) {
            for (count j = 0; j < seq[i]; ++j) {
                urn.push_back(i);
            }
            remaining[i] = seq[i];
            if (seq[i] == 0) {
                --remainingNodes;
            }
        }

        count remainingSize = urn.size();
        index curIndex = 0;

        // only indices that are larger than curIndex will be sampled
        auto sampleAvailableIndex = [&remainingSize, &urn, &alreadyTaken, &curIndex]() -> index {
            index result;
            do {
                --remainingSize;
                index i;
                do {
                    i = Aux::Random::integer(remainingSize);
                    result = urn[i];
                } while (alreadyTaken[result]); // only considering nodes that are not yet connected
                                                // to curIndex
                urn[i] = urn[remainingSize];
            } while (result
                     <= curIndex); // when the sampled nodeID is smaller or equal to the
                                   // curIndex, it is removed from the urn, only if it is greater
                                   // than the curIndex it is available and therefore gets returned
            return result;
        };

        std::vector<index> neighbors;
        for (curIndex = 0; curIndex < seq.size(); ++curIndex) {
            if (remaining[curIndex] == 0)
                continue;
            --remainingNodes;
            if (remaining[curIndex] > remainingNodes)
                break; // Not enough neighbours to complete
            neighbors.clear();
            for (node j = 0; j < remaining[curIndex]; ++j) {
                index availableIndex = sampleAvailableIndex();
                result.addEdge(curIndex, availableIndex);
                neighbors.push_back(availableIndex);
                alreadyTaken[availableIndex] = true;
                --remaining[availableIndex];
                if (remaining[availableIndex] == 0)
                    --remainingNodes;
            }
            for (index a : neighbors) {
                alreadyTaken[a] = false;
            }
        }

        if (curIndex != seq.size())
            continue;
        return result;
    } while (true);
}

} // namespace NetworKit
