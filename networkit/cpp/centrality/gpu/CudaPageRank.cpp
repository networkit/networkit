// networkit-format
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cuda_runtime.h>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

#include "CudaPageRankKernels.hpp"
#include "PageRankBackends.hpp"
#include "PageRankMatrix.hpp"

namespace NetworKit {
namespace GPU {

namespace {

void check(cudaError_t error) {
    if (error != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA error: ") + cudaGetErrorString(error));
    }
}

bool hasDevice() {
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess) {
        // Clears the error, so that it does not show up in later calls.
        cudaGetLastError();
        return false;
    }
    return devices > 0;
}

// Owns an array in device memory.
template <class T>
class DeviceArray {
public:
    explicit DeviceArray(count size) : size(size) {
        check(cudaMalloc(&ptr, std::max<count>(size, 1) * sizeof(T)));
    }

    ~DeviceArray() { cudaFree(ptr); }

    DeviceArray(const DeviceArray &) = delete;
    DeviceArray &operator=(const DeviceArray &) = delete;

    T *data() const { return ptr; }

    void copyFrom(const std::vector<T> &host) {
        check(cudaMemcpy(ptr, host.data(), size * sizeof(T), cudaMemcpyHostToDevice));
    }

    void copyTo(std::vector<T> &host) const {
        host.resize(size);
        check(cudaMemcpy(host.data(), ptr, size * sizeof(T), cudaMemcpyDeviceToHost));
    }

private:
    T *ptr = nullptr;
    count size;
};

} // namespace

std::optional<PageRankRun> cudaPageRank(const Graph &G, const std::vector<uint32_t> &rows,
                                        double damp, double tol, count maxIterations) {
    Aux::SignalHandler handler;
    const count n = G.numberOfNodes();

    if (!hasDevice()) {
        return std::nullopt;
    }

    std::vector<uint64_t> rowBegins(n + 1);
    const count numEntries = fillRowBegins(G, rows, rowBegins.data());

    size_t freeBytes = 0;
    size_t totalBytes = 0;
    check(cudaMemGetInfo(&freeBytes, &totalBytes));
    const count bytes = (n + 1) * sizeof(uint64_t)
                        + numEntries * (sizeof(uint32_t) + sizeof(double)) + 2 * n * sizeof(double)
                        + sizeof(double);
    if (bytes > freeBytes) {
        return std::nullopt;
    }

    std::vector<uint32_t> columns(numEntries);
    std::vector<double> values(numEntries);
    fillEntries(G, damp, rows, rowBegins.data(), columns.data(), values.data());

    DeviceArray<uint64_t> deviceRowBegins(n + 1);
    DeviceArray<uint32_t> deviceColumns(numEntries);
    DeviceArray<double> deviceValues(numEntries);
    std::array<DeviceArray<double>, 2> deviceScores{DeviceArray<double>(n), DeviceArray<double>(n)};
    DeviceArray<double> deviceSquaredChange(1);

    int device = 0;
    cudaDeviceProp properties{};
    check(cudaGetDevice(&device));
    check(cudaGetDeviceProperties(&properties, device));
    INFO("GPUPageRank: running on CUDA, device ", properties.name);

    deviceRowBegins.copyFrom(rowBegins);
    deviceColumns.copyFrom(columns);
    deviceValues.copyFrom(values);
    deviceScores[0].copyFrom(std::vector<double>(n, 1. / static_cast<double>(n)));

    const double teleport = (1. - damp) / static_cast<double>(n);
    count iterations = 0;
    bool isConverged = false;
    do {
        handler.assureRunning();
        check(cudaMemset(deviceSquaredChange.data(), 0, sizeof(double)));
        CudaPageRankKernels::iterate(
            deviceRowBegins.data(), deviceColumns.data(), deviceValues.data(), teleport,
            deviceScores[iterations % 2].data(), deviceScores[(iterations + 1) % 2].data(),
            static_cast<uint32_t>(n), deviceSquaredChange.data());
        check(cudaGetLastError());

        // Waits for the kernel to finish.
        double squaredChange = 0.;
        check(cudaMemcpy(&squaredChange, deviceSquaredChange.data(), sizeof(double),
                         cudaMemcpyDeviceToHost));

        ++iterations;
        isConverged = iterations >= maxIterations || std::sqrt(squaredChange) <= tol;
    } while (!isConverged);

    PageRankRun result;
    deviceScores[iterations % 2].copyTo(result.scores);
    result.iterations = iterations;
    return result;
}

} // namespace GPU
} // namespace NetworKit
