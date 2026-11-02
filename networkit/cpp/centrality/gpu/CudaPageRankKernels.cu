// networkit-format
#include <cstdint>

#include "CudaPageRankKernels.hpp"

namespace NetworKit {
namespace CudaPageRankKernels {

namespace {

constexpr unsigned threadsPerBlock = 256;
constexpr unsigned threadsPerWarp = 32;
constexpr unsigned fullWarp = 0xffffffff;

__device__ double warpSum(double value) {
    for (unsigned offset = threadsPerWarp / 2; offset > 0; offset /= 2) {
        value += __shfl_down_sync(fullWarp, value, offset);
    }
    return value;
}

// One thread per row. The squared changes are summed per block, so that there is one atomic
// addition per block instead of one per row.
__global__ void pageRankIteration(const uint64_t *rowBegins, const uint32_t *columns,
                                  const double *values, double teleport, const double *input,
                                  double *output, uint32_t n, double *squaredChange) {
    __shared__ double warpChanges[threadsPerBlock / threadsPerWarp];

    const uint64_t u = static_cast<uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    double change = 0.;
    if (u < n) {
        double rank = 0.;
        for (auto entry = rowBegins[u]; entry < rowBegins[u + 1]; ++entry) {
            rank += values[entry] * input[columns[entry]];
        }
        rank += teleport;
        output[u] = rank;
        change = (rank - input[u]) * (rank - input[u]);
    }

    const unsigned lane = threadIdx.x % threadsPerWarp;
    const unsigned warp = threadIdx.x / threadsPerWarp;
    change = warpSum(change);
    if (lane == 0) {
        warpChanges[warp] = change;
    }
    __syncthreads();

    if (warp == 0) {
        change = warpSum(lane < threadsPerBlock / threadsPerWarp ? warpChanges[lane] : 0.);
        if (lane == 0) {
            atomicAdd(squaredChange, change);
        }
    }
}

} // namespace

void iterate(const uint64_t *rowBegins, const uint32_t *columns, const double *values,
             double teleport, const double *input, double *output, uint32_t n,
             double *squaredChange) {
    const unsigned blocks = (n + threadsPerBlock - 1) / threadsPerBlock;
    pageRankIteration<<<blocks, threadsPerBlock>>>(rowBegins, columns, values, teleport, input,
                                                   output, n, squaredChange);
}

} // namespace CudaPageRankKernels
} // namespace NetworKit
