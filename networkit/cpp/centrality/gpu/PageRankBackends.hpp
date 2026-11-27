// networkit-format
#ifndef NETWORKIT_CENTRALITY_GPU_PAGE_RANK_BACKENDS_HPP_
#define NETWORKIT_CENTRALITY_GPU_PAGE_RANK_BACKENDS_HPP_

#include <cstdint>
#include <optional>
#include <vector>

#include <networkit/graph/Graph.hpp>

// Each backend is only defined if NetworKit is built with it.
namespace NetworKit {
namespace GPU {

struct PageRankRun {
    // Indexed by the rows of the matrix, see matrixRows.
    std::vector<double> scores;
    count iterations = 0;
};

// Each returns nothing if there is no GPU, or the graph does not fit in its memory. rows are the
// rows of the matrix, see matrixRows.
std::optional<PageRankRun> cudaPageRank(const Graph &G, const std::vector<uint32_t> &rows,
                                        double damp, double tol, count maxIterations);
std::optional<PageRankRun> metalPageRank(const Graph &G, const std::vector<uint32_t> &rows,
                                         double damp, double tol, count maxIterations);

} // namespace GPU
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GPU_PAGE_RANK_BACKENDS_HPP_
