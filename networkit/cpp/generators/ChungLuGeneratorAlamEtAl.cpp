/*
 * ChungLuGeneratorAlamEtAl.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth/Hering
 */

#include <iostream>
#include <numeric>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ChungLuGeneratorAlamEtAl.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

ChungLuGeneratorAlamEtAl::ChungLuGeneratorAlamEtAl(const std::vector<count> &degreeSequence,
                                                   bool parallel)
    : parallel(parallel) {
    n = (count)degreeSequence.size();
    sum_deg = 0;
    // Since the degrees cannot be larger than n, counting sort provides linear runtime complexity
    std::vector<count> buckets{};
    buckets.resize(n, 0);
    for (index i = 0; i < n; i++) {
        count degrees = std::min(degreeSequence[i], n - 1);
        buckets[degrees]++;
        sum_deg += degrees;
    }
    count maxGroupSize = static_cast<count>(
        10 + n / std::sqrt(omp_get_max_threads()) / 2); // Make sure that the parallelization
                                                        // step has enough groups to work with.
    if (!parallel) {
        maxGroupSize = 0xffffffffffffffffL;
    }

    auto fill_groups = [&](index i, count startIndex) {
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
    };

    count startIndex = 0;
    index i = 0;
    while (startIndex < n) {
        while (buckets[i] == 0)
            i++;
        if (parallel && buckets[i] > maxGroupSize) { // this can only happen in the parallel case
                                                     // (for sequential maxGroupSize is ULONG_MAX)
            fill_groups(i, startIndex);
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
    if (p == 0.0)
        return;

    index x = 0 - 1;
    double logP = 1.0 / std::log(1.0 - p);
    bool end_arrived = false;
    do {
        double randVal = real(rng);
        double l = std::floor(std::log(randVal) * logP);
        if (l == std::numeric_limits<double>::infinity())
            continue;
        x += static_cast<index>(l + 1);
        if (x >= end) {
            end_arrived = true;
            continue;
        }

        index u;
        index v;
        if (i == j) {
            u = static_cast<index>(std::floor((1.0 + std::sqrt(1.0 + 8.0 * x)) / 2.0));
            v = x - u * (u - 1) / 2;
        } else {
            u = x / groups[j].vertexCount;
            v = x % groups[j].vertexCount;
        }
        addEdge(groups[i].startIndex + u, groups[j].startIndex + v);
    } while (!end_arrived);
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

    for (index i = 0; i < static_cast<index>(groups.size()); i++) {
        const double p_i = static_cast<double>(groups[i].degrees)
                           * static_cast<double>(groups[i].degrees) / static_cast<double>(sum_deg);
        for (index j = 0; j < i; j++) {
            const double p_j = static_cast<double>(groups[i].degrees)
                               * static_cast<double>(groups[j].degrees)
                               / static_cast<double>(sum_deg);
            edgeSkipping(rng, addEdge, i, j, p_j, groups[i].vertexCount * groups[j].vertexCount);
        }
        edgeSkipping(rng, addEdge, i, i, p_i,
                     groups[i].vertexCount * (groups[i].vertexCount - 1) / 2);
    }

    return G;
}

Graph ChungLuGeneratorAlamEtAl::generateParallel() {
    GraphBuilder builder(n);

    auto addEdge = [&builder](index u, index v) { builder.addHalfEdge(u, v); };

    int max = groups.size() * (groups.size() + 1) / 2;
#pragma omp parallel for schedule(guided)
    for (int xInv = 0; xInv < max; xInv++) {
        int x = max - 1 - xInv; // Starting with the small chunks improves guided scheduling.
        std::mt19937_64 &rng = Aux::Random::getURNG();
        int i = std::floor((-1.0 + std::sqrt(1.0 + 8.0 * x)) / 2.0);
        int j = x - i * (i + 1) / 2;
        if (i != j) {
            const double p_j = static_cast<double>(groups[i].degrees)
                               * static_cast<double>(groups[j].degrees)
                               / static_cast<double>(sum_deg);
            edgeSkipping(rng, addEdge, i, j, p_j, groups[i].vertexCount * groups[j].vertexCount);
        } else {
            const double p_i = static_cast<double>(groups[i].degrees)
                               * static_cast<double>(groups[i].degrees)
                               / static_cast<double>(sum_deg);
            edgeSkipping(rng, addEdge, i, i, p_i,
                         groups[i].vertexCount * (groups[i].vertexCount - 1) / 2);
        }
    }

    return builder.completeGraph();
}

} /* namespace NetworKit */
