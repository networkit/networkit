/*
 * Generator.h
 *
 *  Created on: 03. May 2019
 *      Author: Christopher Weyand <Christopher.Weyand@hpi.de>, Manuel Penschuck <networkit@manuel.jetzt>
 *
 * Code is adopted from https://github.com/chistopher/girgs
 */

#include <functional>
#include <mutex>

#include <omp.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphBuilder.hpp>

#include "Generator.hpp"
#include "SpatialTree.hpp"
#include "WeightScaling.hpp"

namespace NetworKit {
namespace girgs {

std::vector<double> generateWeights(int n, double ple, bool parallel) {
    const auto threads = parallel ? std::max(1, std::min(omp_get_max_threads(), n / 10000)) : 1;
    auto result = std::vector<double>(n);

    #pragma omp parallel num_threads(threads)
    {
        auto& gen = Aux::Random::getURNG();
        auto dist = std::uniform_real_distribution<>{};

        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) {
            result[i] = std::pow((std::pow(0.5*n, -ple + 1) - 1) * dist(gen) + 1, 1 / (-ple + 1));
        }
    }

    return result;
}

std::vector<std::vector<double>> generatePositions(int n, int dimension, bool parallel) {
    const auto threads = parallel ? std::max(1, std::min(omp_get_max_threads(), n / 10000)) : 1;
    auto result = std::vector<std::vector<double>>(n, std::vector<double>(dimension));

    #pragma omp parallel num_threads(threads)
    {
        auto& gen = Aux::Random::getURNG();
        auto dist = std::uniform_real_distribution<>{};

        #pragma omp for schedule(static)
        for(int i=0; i<n; ++i)
            for (int d=0; d<dimension; ++d)
                result[i][d] = dist(gen);
    }

    return result;
}

double scaleWeights(std::vector<double>& weights, double desiredAvgDegree, int dimension, double alpha) {
    // estimate scaling with binary search
    double scaling;
    if(alpha > 8.0)
        scaling = estimateWeightScalingThreshold(weights, desiredAvgDegree, dimension);
    else if(alpha > 0.0 && alpha != 1.0)
        scaling = estimateWeightScaling(weights, desiredAvgDegree, dimension, alpha);
    else
        throw("Scaling requires alpha > 0 with alpha != 1");

    // scale weights
    for(auto& each : weights)
        each *= scaling;
    return scaling;
}

Graph generateEdges(std::vector<double>& weights, std::vector<std::vector<double>>& positions, double alpha, bool keep_input) {

    assert(weights.size() == positions.size());

    // undirected, unweighted graph with n nodes
    GraphBuilder builder(weights.size(), false, false);

    using edge_vector = std::vector<std::pair<int, int>>;
    edge_vector result;

    std::vector<std::pair<
        edge_vector,
        uint64_t[31] /* avoid false sharing */
    > > local_edges(omp_get_max_threads());

    constexpr auto block_size = size_t{1} << 20;

    std::mutex m;
    auto flush = [&] (const edge_vector& local) {
        std::lock_guard<std::mutex> lock(m);
        for(auto edge : local)
            builder.addHalfEdge(edge.first, edge.second); // GraphBuilder's interface is not really threadsafe :(
    };

    auto addEdge = [&](int u, int v, int tid) {
        auto& local = local_edges[tid].first;
        local.emplace_back(u,v);
        if (local.size() == block_size) {
            flush(local);
            local.clear();
        }
    };
    auto dimension = positions.front().size();

    auto maybe_clear_input = [&] {
        if (keep_input)
            return;

        weights.clear();
        positions.clear();
    };

    switch(dimension) {
        case 1: {auto tree = makeSpatialTree<1>(weights, positions, alpha, addEdge); maybe_clear_input(); tree.generateEdges(); break;}
        case 2: {auto tree = makeSpatialTree<2>(weights, positions, alpha, addEdge); maybe_clear_input(); tree.generateEdges(); break;}
        case 3: {auto tree = makeSpatialTree<3>(weights, positions, alpha, addEdge); maybe_clear_input(); tree.generateEdges(); break;}
        case 4: {auto tree = makeSpatialTree<4>(weights, positions, alpha, addEdge); maybe_clear_input(); tree.generateEdges(); break;}
        case 5: {auto tree = makeSpatialTree<5>(weights, positions, alpha, addEdge); maybe_clear_input(); tree.generateEdges(); break;}
        default:
            throw std::runtime_error("Support only dimensions 1 to 5");
    }

    for(const auto& v : local_edges)
        flush(v.first);


    return builder.toGraph(true);
}

} // namespace girgs
} // namespace NetworKit

