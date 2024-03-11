/*
 * ChungLu.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include <iostream>
#include <numeric>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ChungLuGeneratorAlamEtAl.hpp>
#include <networkit/graph/ParallelGraphBuilder.hpp>

namespace NetworKit {

ChungLuGeneratorAlamEtAl::ChungLuGeneratorAlamEtAl(const std::vector<count> &degreeSequence,
                                                   bool parallel)
    : parallel(parallel) {
    n = (count)degreeSequence.size();
    sum_deg = 0;
    // Since the degrees cannot be larger than n, a bucket sort is O(n) here:
    std::vector<count> buckets{};
    buckets.resize(n, 0);
    for (index i = 0; i < n; i++) {
        count degrees = std::min(degreeSequence[i], n - 1);
        buckets[degrees]++;
        sum_deg += degrees;
    }
    count maxGroupSize =
        (count)(10
                + n / std::sqrt(omp_get_max_threads()) / 2); // Make sure that the parallelization
                                                             // step has enough groups to work with.
    if (!parallel) {
        maxGroupSize = 0xffffffffffffffffL;
    }
    count startIndex = 0;
    index i = 0;
    while (startIndex < n) {
        while (buckets[i] == 0)
            i++;
        if (buckets[i] > maxGroupSize) {
            count pieces = buckets[i] / maxGroupSize + 1;
            count pieceSize = buckets[i] / pieces;
            count leftOver = buckets[i] % pieces;
            for (index piece = 0; piece < pieces; piece++) {
                count size = pieceSize;
                if (piece < leftOver) {
                    size++;
                }
                groups.push_back(ChungLuGeneratorAlamEtAl::VertexGroup{i, size, startIndex});
                startIndex += size;
            }
        } else {
            groups.push_back(ChungLuGeneratorAlamEtAl::VertexGroup{i, buckets[i], startIndex});
            startIndex += buckets[i];
        }
        i++;
    }
    groups.shrink_to_fit();
}

double real(std::mt19937_64 &rng) {
    thread_local static std::uniform_real_distribution<double> dist{};
    return dist(rng);
}

template <typename F>
void ChungLuGeneratorAlamEtAl::edgeSkipping(std::mt19937_64 &rng, F &&addEdge, index i, index j,
                                            double p, index end) {
    index x = 0 - 1;
    double logP = 1.0 / std::log(1.0 - p);
    while (true) {
        double randVal = real(rng);
        double l = std::floor(std::log(randVal) * logP);
        x += (index)l + 1;
        if (x >= end)
            break;

        index u;
        index v;
        if (i == j) {
            u = (index)std::floor((1.0 + std::sqrt(1.0 + 8.0 * x)) / 2.0);
            v = x - u * (u - 1) / 2;
        } else {
            u = x / groups[j].vertexCount;
            v = x % groups[j].vertexCount;
        }
        addEdge(groups[i].startIndex + u, groups[j].startIndex + v);
    }
}

Graph ChungLuGeneratorAlamEtAl::generate() {
    if (parallel) {
        return generateParallel();
    } else {
        return generateSequential();
    }
}

Graph ChungLuGeneratorAlamEtAl::generateSequential() {
    Graph G(n);
    std::mt19937_64 &rng = Aux::Random::getURNG();

    auto addEdge = [&G](index u, index v) { G.addEdge(u, v); };

    for (index i = 0; i < (index)groups.size(); i++) {
        for (index j = 0; j < i; j++) {
            edgeSkipping(rng, addEdge, i, j,
                         (double)groups[i].degrees * (double)groups[j].degrees / (double)sum_deg,
                         groups[i].vertexCount * groups[j].vertexCount);
        }
        edgeSkipping(rng, addEdge, i, i,
                     (double)groups[i].degrees * (double)groups[i].degrees / (double)sum_deg,
                     groups[i].vertexCount * (groups[i].vertexCount - 1) / 2);
    }

    return G;
}

Graph ChungLuGeneratorAlamEtAl::generateParallel() {
    ParallelGraphBuilder builder(n);

    auto addEdge = [&builder](index u, index v) { builder.addEdgeParallel(u, v); };

    int max = groups.size() * (groups.size() + 1) / 2;
#pragma omp parallel for schedule(guided)
    for (int xInv = 0; xInv < max; xInv++) {
        int x = max - 1 - xInv; // Starting with the small chunks improves guided scheduling.
        std::mt19937_64 &rng = Aux::Random::getURNG();
        int i = std::floor((-1.0 + std::sqrt(1.0 + 8.0 * x)) / 2.0);
        int j = x - i * (i + 1) / 2;
        if (i != j) {
            edgeSkipping(rng, addEdge, i, j,
                         (double)groups[i].degrees * (double)groups[j].degrees / (double)sum_deg,
                         groups[i].vertexCount * groups[j].vertexCount);
        } else {
            edgeSkipping(rng, addEdge, i, i,
                         (double)groups[i].degrees * (double)groups[i].degrees / (double)sum_deg,
                         groups[i].vertexCount * (groups[i].vertexCount - 1) / 2);
        }
    }

    return builder.toGraphParallel();
}

} /* namespace NetworKit */