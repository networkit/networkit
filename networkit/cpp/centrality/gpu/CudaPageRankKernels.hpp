// networkit-format
#ifndef NETWORKIT_CENTRALITY_CUDA_PAGE_RANK_KERNELS_HPP_
#define NETWORKIT_CENTRALITY_CUDA_PAGE_RANK_KERNELS_HPP_

#include <cstdint>

namespace NetworKit {
namespace CudaPageRankKernels {

/**
 * Launches one PageRank iteration on the n rows of the matrix: output = teleport + matrix * input.
 * Adds the squared L2 norm of output - input to squaredChange. All pointers are in device memory.
 */
void iterate(const uint64_t *rowBegins, const uint32_t *columns, const double *values,
             double teleport, const double *input, double *output, uint32_t n,
             double *squaredChange);

} // namespace CudaPageRankKernels
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_CUDA_PAGE_RANK_KERNELS_HPP_
