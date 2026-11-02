// networkit-format
#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/centrality/GPUPageRank.hpp>
#include <networkit/centrality/PageRank.hpp>

#include "gpu/PageRankBackends.hpp"
#include "gpu/PageRankMatrix.hpp"

namespace NetworKit {

namespace {

// The build defines NETWORKIT_WITH_CUDA and NETWORKIT_WITH_METAL for this file only.
std::optional<GPU::PageRankRun> runBackend(GPUBackend backend, [[maybe_unused]] const Graph &G,
                                           [[maybe_unused]] const std::vector<uint32_t> &rows,
                                           [[maybe_unused]] double damp,
                                           [[maybe_unused]] double tol,
                                           [[maybe_unused]] count maxIterations) {
    switch (backend) {
    case GPUBackend::CUDA:
#ifdef NETWORKIT_WITH_CUDA
        return GPU::cudaPageRank(G, rows, damp, tol, maxIterations);
#else
        return std::nullopt;
#endif
    case GPUBackend::METAL:
#ifdef NETWORKIT_WITH_METAL
        return GPU::metalPageRank(G, rows, damp, tol, maxIterations);
#else
        return std::nullopt;
#endif
    default:
        return std::nullopt;
    }
}

} // namespace

GPUPageRank::GPUPageRank(const Graph &G, double damp, double tol, GPUBackend backend)
    : Centrality(G, true), damp(damp), tol(tol), backend(backend) {
    if (!isBuilt(backend)) {
        throw std::runtime_error("NetworKit is built without this GPU backend.");
    }
}

bool GPUPageRank::isBuilt(GPUBackend backend) {
    switch (backend) {
    case GPUBackend::CUDA:
#ifdef NETWORKIT_WITH_CUDA
        return true;
#else
        return false;
#endif
    case GPUBackend::METAL:
#ifdef NETWORKIT_WITH_METAL
        return true;
#else
        return false;
#endif
    default:
        return true;
    }
}

void GPUPageRank::run() {
    const count n = G.numberOfNodes();

    scoreData.assign(G.upperNodeIdBound(), 0.);
    max = 0.;
    iterations = 0;
    used = GPUBackend::CPU;

    if (n == 0) {
        hasRun = true;
        return;
    }

    // The backends index rows and columns with 32 bits.
    std::optional<GPU::PageRankRun> result;
    std::vector<uint32_t> rows;
    if (backend != GPUBackend::CPU && n < std::numeric_limits<uint32_t>::max()) {
        rows = GPU::matrixRows(G);
        for (const auto candidate : {GPUBackend::CUDA, GPUBackend::METAL}) {
            if ((backend == GPUBackend::AUTO || backend == candidate) && isBuilt(candidate)) {
                result = runBackend(candidate, G, rows, damp, tol, maxIterations);
                if (result) {
                    used = candidate;
                    break;
                }
            }
        }
    }

    if (!result) {
        if (backend == GPUBackend::CPU) {
            INFO("GPUPageRank: running on the CPU");
        } else {
            INFO("GPUPageRank: no usable GPU, running on the CPU");
        }
        runOnCpu();
        return;
    }

    iterations = result->iterations;
    G.parallelForNodes([&](node u) -> void { scoreData[u] = result->scores[rows[u]]; });

    // Without sink handling, the scores of a directed graph do not sum to 1.
    if (G.isDirected()) {
        const double sum = G.parallelSumForNodes([&](node u) -> double { return scoreData[u]; });
        G.parallelForNodes([&](node u) -> void { scoreData[u] /= sum; });
    }

    max = *std::max_element(scoreData.begin(), scoreData.end());
    hasRun = true;
}

void GPUPageRank::runOnCpu() {
    PageRank pageRank(G, damp, tol);
    pageRank.maxIterations = maxIterations;
    pageRank.run();

    scoreData = pageRank.scores();
    max = pageRank.maximum();
    iterations = pageRank.numberOfIterations();
    hasRun = true;
}

double GPUPageRank::maximum() {
    return max;
}

} /* namespace NetworKit */
