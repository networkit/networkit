
#include <networkit/generators/ConfigurationModel.hpp>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/graph/Graph.hpp>

#include <iostream>

namespace NetworKit {

ConfigurationModel::ConfigurationModel(const std::vector<count> &sequence)
    : StaticDegreeSequenceGenerator(sequence) {

    // Throw an error if the sequence isn't realizable:
    std::sort(seq.begin(), seq.end(), std::greater<count>());
    HavelHakimiGenerator(seq).generate();
}

Graph ConfigurationModel::generate() {
    do {
        Graph result{seq.size()};

        std::vector<index> urn{};
        std::vector<count> remaining(seq.size());
        std::vector<bool> alreadyTaken(seq.size(), false);
        index remainingNodes = seq.size();

        // Fill the urn:
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

        auto sample = [&remainingSize, &urn, &alreadyTaken, &curIndex]() -> index {
            index result;
            do {
                --remainingSize;
                index i;
                do {
                    i = Aux::Random::integer(remainingSize);
                    result = urn[i];
                } while (alreadyTaken[result]);
                urn[i] = urn[remainingSize];
            } while (result <= curIndex);
            return result;
        };

        std::vector<index> neighbors;
        for (curIndex = 0; curIndex < seq.size(); ++curIndex) {
            if (remaining[curIndex] == 0)
                continue;
            remainingNodes--;
            if (remaining[curIndex] > remainingNodes)
                break; // Impossible
            neighbors.clear();
            for (node j = 0; j < remaining[curIndex]; ++j) {
                index a = sample();
                result.addEdge(curIndex, a);
                neighbors.push_back(a);
                alreadyTaken[a] = true;
                --remaining[a];
                if (remaining[a] == 0)
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
