
#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include <networkit/graph/Graph.hpp>
#include <networkit/scd/ApproximatePageRank.hpp>

namespace NetworKit {
TEST(ApproximatePageRankGTest, RearmsAtRoundedEpsilonBoundary) {
    Graph graph(34, false, false);
    const std::pair<node, node> edges[] = {
        {0, 1},   {0, 2},   {0, 3},   {0, 4},   {0, 5},   {0, 6},   {0, 7},   {0, 8},   {0, 10},
        {0, 11},  {0, 12},  {0, 13},  {0, 17},  {0, 19},  {0, 21},  {0, 31},  {1, 2},   {1, 3},
        {1, 7},   {1, 13},  {1, 17},  {1, 19},  {1, 21},  {1, 30},  {2, 3},   {2, 7},   {2, 8},
        {2, 9},   {2, 13},  {2, 27},  {2, 28},  {2, 32},  {3, 7},   {3, 12},  {3, 13},  {4, 6},
        {4, 10},  {5, 6},   {5, 10},  {5, 16},  {6, 16},  {8, 30},  {8, 32},  {8, 33},  {9, 33},
        {13, 33}, {14, 32}, {14, 33}, {15, 32}, {15, 33}, {18, 32}, {18, 33}, {19, 33}, {20, 32},
        {20, 33}, {22, 32}, {22, 33}, {23, 25}, {23, 27}, {23, 29}, {23, 32}, {23, 33}, {24, 25},
        {24, 27}, {24, 31}, {25, 31}, {26, 29}, {26, 33}, {27, 33}, {28, 31}, {28, 33}, {29, 32},
        {29, 33}, {30, 32}, {30, 33}, {31, 32}, {31, 33}, {32, 33}};
    for (auto [u, v] : edges)
        graph.addEdge(u, v);
    // Solve the lazy-walk equations independently by power iteration in
    // extended precision, well below the tested approximation tolerance.
    const long double alpha = static_cast<long double>(0x1.3333333333333p-3);
    std::vector<long double> reference(34, 0.0L);
    for (unsigned step = 0; step < 2000; ++step) {
        std::vector<long double> next(34, 0.0L);
        next[0] = alpha;
        for (node u = 0; u < 34; ++u) {
            auto mass = (1 - alpha) * reference[u] / 2;
            next[u] += mass;
            graph.forNeighborsOf(u, [&](node v) { next[v] += mass / graph.degree(u); });
        }
        reference.swap(next);
    }
    const double boundary = 0x1.dc5b7a26f5e67p-10;
    for (double eps : {std::nextafter(boundary, 0.0), boundary, std::nextafter(boundary, 1.0)}) {
        ApproximatePageRank algorithm(graph, static_cast<double>(alpha), eps);
        std::vector<long double> actual(34, 0.0L);
        for (auto [u, value] : algorithm.run(0))
            actual[u] = value;
        for (node u = 0; u < 34; ++u) {
            EXPECT_GE(actual[u], reference[u] - eps * graph.degree(u) - 1e-12L) << u;
            EXPECT_LE(actual[u], reference[u] + 1e-12L) << u;
        }
    }
}
} // namespace NetworKit
