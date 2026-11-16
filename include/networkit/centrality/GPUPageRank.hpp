// networkit-format
#ifndef NETWORKIT_CENTRALITY_GPU_PAGE_RANK_HPP_
#define NETWORKIT_CENTRALITY_GPU_PAGE_RANK_HPP_

#include <limits>

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

enum class GPUBackend {
    // The first backend that finds a GPU, in the order CUDA, Metal, then the CPU.
    AUTO,
    // NVIDIA GPUs, in double precision. Needs NETWORKIT_WITH_CUDA=ON.
    CUDA,
    // Apple GPUs, in single precision. Needs NETWORKIT_WITH_METAL=ON, on macOS only.
    METAL,
    // PageRank, on the CPU.
    CPU
};

/**
 * @ingroup centrality
 * Prototype of PageRank that runs on a GPU, with the default parameters of PageRank: no sink
 * handling, no normalization, and the L2 norm as stopping criterion. Which GPUs it supports depends
 * on the backends NetworKit is built with, see GPUBackend.
 *
 * Metal has no double, so its scores are only as accurate as single precision allows, whatever the
 * tolerance. Its iteration also stops once the change of the scores is down to rounding errors.
 *
 * If the backend finds no GPU, or the graph does not fit in its memory, this runs PageRank on the
 * CPU.
 */
class GPUPageRank final : public Centrality {

public:
    /**
     * Constructs the GPUPageRank class for the Graph @a G
     *
     * @param[in] G Graph to be processed.
     * @param[in] damp Damping factor of the PageRank algorithm.
     * @param[in] tol Error tolerance for PageRank iteration.
     * @param[in] backend Backend to run on. Throws if NetworKit is built without it.
     */
    GPUPageRank(const Graph &G, double damp = 0.85, double tol = 1e-8,
                GPUBackend backend = GPUBackend::AUTO);

    /**
     * Computes page rank on the graph passed in constructor.
     */
    void run() override;

    /**
     * Returns the maximum PageRank score
     */
    double maximum() override;

    /**
     * Return the number of iterations performed by the algorithm.
     *
     * @return Number of iterations performed by the algorithm.
     */
    count numberOfIterations() const {
        assureFinished();
        return iterations;
    }

    /**
     * Return the backend the algorithm ran on: CPU if it fell back to the CPU, never AUTO.
     *
     * @return The backend the algorithm ran on.
     */
    GPUBackend usedBackend() const {
        assureFinished();
        return used;
    }

    /**
     * Return whether NetworKit is built with the given backend. AUTO and CPU are always built.
     *
     * @param[in] backend The backend.
     * @return True if NetworKit is built with @a backend.
     */
    static bool isBuilt(GPUBackend backend);

    // Maximum number of iterations allowed
    count maxIterations = std::numeric_limits<count>::max();

private:
    void runOnCpu();

    double damp;
    double tol;
    GPUBackend backend;
    count iterations = 0;
    GPUBackend used = GPUBackend::CPU;
    double max = 0.;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_GPU_PAGE_RANK_HPP_
