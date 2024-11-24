/*
 * RmatGenerator.cpp
 *
 *  Created on: 18.03.2014
 *      Author: Henning, cls
 *
 * Uses the algorithm described by Hübschle-Schneider and Sanders in
 * "Linear Work Generation of R-MAT Graphs" https://arxiv.org/abs/1905.03525
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/RmatGenerator.hpp>

namespace NetworKit {

struct Entry {
    uint32_t i;
    uint32_t j;
    uint32_t numberOfBits;
    double probability;
};

bool priorityQueueCompare(const Entry &a, const Entry &b) {
    if (a.numberOfBits == 32 || b.numberOfBits == 32) {
        return a.numberOfBits > b.numberOfBits;
    }
    return a.probability < b.probability;
}

// std::priority_queue doesn't support iterating over its elements.
struct PriorityQueue {
    std::vector<Entry> &internal;

    PriorityQueue(std::vector<Entry> &vec) : internal(vec) {}

    void push(Entry &&entry) {
        internal.push_back(entry);
        std::push_heap(internal.begin(), internal.end(), priorityQueueCompare);
    }

    Entry pop() {
        std::pop_heap(internal.begin(), internal.end(), priorityQueueCompare);
        Entry result = internal.back();
        internal.pop_back();
        return result;
    }
};

void generateEntries(std::vector<Entry> &entryList, count size, double a, double b, double c,
                     double d) {
    entryList.reserve(size);
    PriorityQueue priorityQueue(entryList);
    priorityQueue.push(Entry{0, 0, 1, a});
    priorityQueue.push(Entry{0, 1, 1, b});
    priorityQueue.push(Entry{1, 0, 1, c});
    priorityQueue.push(Entry{1, 1, 1, d});
    while (entryList.size() <= size - 3) {
        // Take the entry with the highest probability and split it up.
        // This increases the average number of bits compared to a uniform distribution.
        Entry old = priorityQueue.pop();
        priorityQueue.push(
            Entry{old.i << 1 | 0, old.j << 1 | 0, old.numberOfBits + 1, a * old.probability});
        priorityQueue.push(
            Entry{old.i << 1 | 0, old.j << 1 | 1, old.numberOfBits + 1, b * old.probability});
        priorityQueue.push(
            Entry{old.i << 1 | 1, old.j << 1 | 0, old.numberOfBits + 1, c * old.probability});
        priorityQueue.push(
            Entry{old.i << 1 | 1, old.j << 1 | 1, old.numberOfBits + 1, d * old.probability});
    }
}

RmatGenerator::RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d,
                             bool weighted, count reduceNodes, bool discardSelfLoops)
    : scale(scale), edgeFactor(edgeFactor), weighted(weighted), reduceNodes(reduceNodes),
      discardSelfLoops(discardSelfLoops) {
    if (scale > 63)
        throw std::runtime_error("Cannot generate more than 2^63 nodes");
    double sum = a + b + c + d;
    INFO("sum of probabilities: ", sum);
    if (!Aux::NumericTools::equal(sum, 1.0, 0.0001))
        throw std::runtime_error("Probabilities in Rmat have to sum to 1.");
    defaultEdgeWeight = 1.0;

    count size = 1 << std::min(static_cast<count>(12),
                               static_cast<count>(std::round(scale * 5.0 / 8.0) + 1));

    assert((size & (size - 1)) == 0); // The size must be a power of two.
    std::vector<Entry> entryList;
    generateEntries(entryList, size, a, b, c, d);
    while (entryList.size() < size) { // Fill the remaining entries with 0-probability entries.
        entryList.push_back(Entry{0, 0, 0, 0});
    }

    // Construct the alias table:
    mask = size - 1;
    bits.resize(size);
    numberOfBits.resize(size);
    coinFlipProbability.resize(size);
    coinFlipReplacement.resize(size);
    for (count i = 0; i < size; i++) {
        bits[i] = std::make_pair(entryList[i].i, entryList[i].j);
        numberOfBits[i] = entryList[i].numberOfBits;
        coinFlipProbability[i] = 0;
        coinFlipReplacement[i] = 0;
    }
    double baseProbability = 1.0 / static_cast<double>(size);
    count lastOverfullIndex = 0;
    count lastUnderfullIndex = 0;
    while (lastUnderfullIndex != size && lastOverfullIndex != size) {
        while (entryList[lastUnderfullIndex].probability >= baseProbability) {
            lastUnderfullIndex++;
            if (lastUnderfullIndex == size)
                return;
        }
        int curUnderfullIndex = lastUnderfullIndex;
        while (lastOverfullIndex != size) {
            while (entryList[lastOverfullIndex].probability <= baseProbability) {
                lastOverfullIndex++;
                if (lastOverfullIndex == size)
                    return;
            }
            double delta = baseProbability - entryList[curUnderfullIndex].probability;
            entryList[curUnderfullIndex].probability = baseProbability;
            entryList[lastOverfullIndex].probability -= delta;
            coinFlipReplacement[curUnderfullIndex] = lastOverfullIndex;
            coinFlipProbability[curUnderfullIndex] = static_cast<uint32_t>(
                delta / baseProbability * std::numeric_limits<uint32_t>::max());
            if (entryList[lastOverfullIndex].probability < baseProbability
                && lastOverfullIndex < lastUnderfullIndex) {
                curUnderfullIndex = lastOverfullIndex;
            } else
                break;
        }
    }
}

std::pair<node, node> RmatGenerator::sampleEdge(uint8_t input_bits) {
    std::pair<node, node> result{0, 0};

    auto sample = [this]() {
        uint64_t randomNumber = Aux::Random::integer();

        uint32_t index = randomNumber & mask;
        uint32_t coinFlip = randomNumber >> 32;
        if (coinFlip <= coinFlipProbability[index]) {
            index = coinFlipReplacement[index];
        }
        curBits = bits[index];
        remainingBits = numberOfBits[index];
    };

    do {
        if (remainingBits >= input_bits) {
            remainingBits -= input_bits;
            result.first = result.first << input_bits | curBits.first >> remainingBits;
            result.second = result.second << input_bits | curBits.second >> remainingBits;
            curBits.first &= (1 << remainingBits) - 1;
            curBits.second &= (1 << remainingBits) - 1;
            return result;
        }
        result.first = result.first << remainingBits | curBits.first;
        result.second = result.second << remainingBits | curBits.second;
        input_bits -= remainingBits;
        sample();
    } while (true);
}

Graph RmatGenerator::generate() {
    double n = std::pow(2.0, static_cast<double>(scale));
    if (n <= reduceNodes) {
        throw std::runtime_error("Error, shall delete more nodes than the graph originally has");
    }
    // when nodes are deleted, all nodes have less neighbors
    count numEdges = n * edgeFactor * n * 1.0 / static_cast<double>(n - reduceNodes);
    count wantedEdges = (n - reduceNodes) * edgeFactor;
    Graph G(n - reduceNodes, weighted);
    // Reset the internal state of the alias table:
    curBits = {0, 0};
    remainingBits = 0;

    Aux::Random::setSeed(Aux::Random::integer(), false);

    if (reduceNodes > 0) {
        std::vector<node> nodemap(n, 0);

        for (count deletedNodes = 0; deletedNodes < reduceNodes;) {
            node u = Aux::Random::index(n);
            if (nodemap[u] == 0) {
                nodemap[u] = none;
                ++deletedNodes;
            }
        }

        for (node i = 0, u = 0; i < n; ++i) {
            if (nodemap[i] == 0) {
                nodemap[i] = u;
                ++u;
            }
        }

        node u, v;
        if (weighted) {
            for (index e = 0; e < numEdges; ++e) {
                std::tie(u, v) = sampleEdge(static_cast<uint8_t>(scale));
                u = nodemap[u];
                v = nodemap[v];
                if (discardSelfLoops && u == v)
                    continue;
                if (u != none && v != none) {
                    G.increaseWeight(u, v, defaultEdgeWeight);
                }
            }
        } else {
            while (G.numberOfEdges() < wantedEdges) {
                std::tie(u, v) = sampleEdge(static_cast<uint8_t>(scale));
                u = nodemap[u];
                v = nodemap[v];
                if (discardSelfLoops && u == v)
                    continue;
                if (u != none && v != none && !G.hasEdge(u, v)) {
                    G.addEdge(u, v);
                }
            }
        }
    } else {
        if (weighted) {
            for (index e = 0; e < numEdges; ++e) {
                std::pair<node, node> drawnEdge = sampleEdge(static_cast<uint8_t>(scale));
                if (discardSelfLoops && drawnEdge.first == drawnEdge.second)
                    continue;
                G.increaseWeight(drawnEdge.first, drawnEdge.second, defaultEdgeWeight);
            }
        } else {
            while (G.numberOfEdges() < wantedEdges) {
                std::pair<node, node> drawnEdge = sampleEdge(static_cast<uint8_t>(scale));
                if (discardSelfLoops && drawnEdge.first == drawnEdge.second)
                    continue;
                if (!G.hasEdge(drawnEdge.first, drawnEdge.second) && G.hasNode(drawnEdge.first)
                    && G.hasNode(drawnEdge.second)) {
                    G.addEdge(drawnEdge.first, drawnEdge.second);
                }
            }
        }
    }

    G.shrinkToFit();
    return G;
}

} /* namespace NetworKit */
