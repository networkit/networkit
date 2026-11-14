#include <metal_stdlib>
using namespace metal;

// One PageRank iteration on row u of the matrix: output = teleport + matrix * input.
kernel void pageRankIteration(const device uint *rowBegins [[buffer(0)]],
                              const device uint *columns [[buffer(1)]],
                              const device float *values [[buffer(2)]],
                              const device float *input [[buffer(3)]],
                              device float *output [[buffer(4)]],
                              constant float &teleport [[buffer(5)]],
                              uint u [[thread_position_in_grid]]) {
    float rank = teleport;
    for (uint entry = rowBegins[u]; entry < rowBegins[u + 1]; ++entry) {
        rank += values[entry] * input[columns[entry]];
    }
    output[u] = rank;
}
