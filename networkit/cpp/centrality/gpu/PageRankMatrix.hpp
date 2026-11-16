// networkit-format
#ifndef NETWORKIT_CENTRALITY_GPU_PAGE_RANK_MATRIX_HPP_
#define NETWORKIT_CENTRALITY_GPU_PAGE_RANK_MATRIX_HPP_

#include <cstdint>
#include <vector>

#include <networkit/graph/Graph.hpp>

// Builds the matrix of the PageRank iteration in CSR format, into memory provided by the backends:
// row u holds damp * weight(v, u) / weightedDegree(v) for each in-neighbor v of u.
namespace NetworKit {
namespace GPU {

// Numbers the rows of the matrix. Removed nodes get no row, so that kernels do not skip them.
inline std::vector<uint32_t> matrixRows(const Graph &G) {
    std::vector<uint32_t> rows(G.upperNodeIdBound(), 0);
    uint32_t nextRow = 0;
    G.forNodes([&](node u) -> void {
        rows[u] = nextRow;
        ++nextRow;
    });
    return rows;
}

// Writes the numberOfNodes() + 1 row offsets of the matrix, and returns its number of entries.
template <class Offset>
count fillRowBegins(const Graph &G, const std::vector<uint32_t> &rows, Offset *rowBegins) {
    rowBegins[0] = 0;
    G.forNodes([&](node u) -> void {
        rowBegins[rows[u] + 1] = rowBegins[rows[u]] + static_cast<Offset>(G.degreeIn(u));
    });
    return static_cast<count>(rowBegins[G.numberOfNodes()]);
}

template <class Offset, class Value>
void fillEntries(const Graph &G, double damp, const std::vector<uint32_t> &rows,
                 const Offset *rowBegins, uint32_t *columns, Value *values) {
    std::vector<double> degrees(G.upperNodeIdBound(), 0.);
    G.parallelForNodes(
        [&](node u) -> void { degrees[u] = static_cast<double>(G.weightedDegree(u)); });
    G.balancedParallelForNodes([&](node u) -> void {
        auto entry = rowBegins[rows[u]];
        G.forInEdgesOf(u, [&](node, node v, edgeweight weight) -> void {
            columns[entry] = rows[v];
            values[entry] = static_cast<Value>(damp * weight / degrees[v]);
            ++entry;
        });
    });
}

} // namespace GPU
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_GPU_PAGE_RANK_MATRIX_HPP_
