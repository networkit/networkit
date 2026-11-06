// networkit-format
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/GPUPageRank.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/RmatGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

namespace {

const std::vector<GPUBackend> allBackends{GPUBackend::AUTO, GPUBackend::CUDA, GPUBackend::METAL,
                                          GPUBackend::CPU};

std::vector<GPUBackend> builtBackends() {
    std::vector<GPUBackend> backends;
    std::copy_if(allBackends.begin(), allBackends.end(), std::back_inserter(backends),
                 [](GPUBackend backend) -> bool { return GPUPageRank::isBuilt(backend); });
    return backends;
}

std::string backendName(GPUBackend backend) {
    switch (backend) {
    case GPUBackend::AUTO:
        return "AUTO";
    case GPUBackend::CUDA:
        return "CUDA";
    case GPUBackend::METAL:
        return "METAL";
    default:
        return "CPU";
    }
}

void expectMatchesPageRank(const Graph &G) {
    PageRank expected(G);
    expected.run();

    for (const auto backend : builtBackends()) {
        SCOPED_TRACE(backendName(backend));
        GPUPageRank gpu(G, 0.85, 1e-8, backend);
        gpu.run();

        if (gpu.usedBackend() == GPUBackend::METAL) {
            // Stopping at a change of 1e-8 leaves both algorithms up to damp / (1 - damp) * 1e-8
            // from the fixed point, and single precision adds an error proportional to the score.
            const auto tolerance = [](double score) -> double { return 1e-7 + 1e-5 * score; };
            EXPECT_NEAR(gpu.maximum(), expected.maximum(), tolerance(expected.maximum()));
            G.forNodes([&](node u) -> void {
                EXPECT_NEAR(gpu.score(u), expected.score(u), tolerance(expected.score(u)));
            });
        } else {
            EXPECT_EQ(gpu.numberOfIterations(), expected.numberOfIterations());
            EXPECT_NEAR(gpu.maximum(), expected.maximum(), 1e-10);
            G.forNodes(
                [&](node u) -> void { EXPECT_NEAR(gpu.score(u), expected.score(u), 1e-10); });
        }
    }
}

} // namespace

class GPUPageRankGTest : public testing::Test {};

TEST_F(GPUPageRankGTest, testUndirected) {
    Aux::Random::setSeed(42, false);
    expectMatchesPageRank(ErdosRenyiGenerator(500, 0.02).generate());
    expectMatchesPageRank(RmatGenerator(9, 8, 0.57, 0.19, 0.19, 0.05).generate());
}

TEST_F(GPUPageRankGTest, testDirected) {
    Aux::Random::setSeed(42, false);
    expectMatchesPageRank(ErdosRenyiGenerator(500, 0.02, true).generate());
}

TEST_F(GPUPageRankGTest, testWeighted) {
    Aux::Random::setSeed(42, false);
    for (const bool directed : {false, true}) {
        auto G = GraphTools::toWeighted(ErdosRenyiGenerator(300, 0.03, directed).generate());
        G.forEdges([&](node u, node v) -> void { G.setWeight(u, v, Aux::Random::real(0.5, 3.)); });
        expectMatchesPageRank(G);
    }
}

TEST_F(GPUPageRankGTest, testRemovedNodes) {
    Aux::Random::setSeed(42, false);
    auto G = ErdosRenyiGenerator(200, 0.05).generate();
    G.removeNode(0);
    G.removeNode(17);
    expectMatchesPageRank(G);
}

TEST_F(GPUPageRankGTest, testSmallGraph) {
    // With Metal, the norm of the scores of so few nodes makes a tolerance of 1e-8 hard to reach
    // with floats, so the iteration has to stop at the rounding errors instead.
    Graph G(4);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    expectMatchesPageRank(G);
}

TEST_F(GPUPageRankGTest, testMaxIterations) {
    Aux::Random::setSeed(42, false);
    const auto G = ErdosRenyiGenerator(100, 0.1).generate();
    for (const auto backend : builtBackends()) {
        SCOPED_TRACE(backendName(backend));
        GPUPageRank gpu(G, 0.85, 1e-8, backend);
        gpu.maxIterations = 3;
        gpu.run();
        EXPECT_EQ(gpu.numberOfIterations(), 3);
    }
}

TEST_F(GPUPageRankGTest, testEmptyGraph) {
    const Graph G;
    for (const auto backend : builtBackends()) {
        SCOPED_TRACE(backendName(backend));
        GPUPageRank gpu(G, 0.85, 1e-8, backend);
        gpu.run();
        EXPECT_EQ(gpu.numberOfIterations(), 0);
    }
}

TEST_F(GPUPageRankGTest, testBackends) {
    const Graph G(10);
    for (const auto backend : allBackends) {
        SCOPED_TRACE(backendName(backend));
        if (!GPUPageRank::isBuilt(backend)) {
            EXPECT_THROW(GPUPageRank(G, 0.85, 1e-8, backend), std::runtime_error);
            continue;
        }

        GPUPageRank gpu(G, 0.85, 1e-8, backend);
        gpu.run();
        EXPECT_NE(gpu.usedBackend(), GPUBackend::AUTO);
        // Without a GPU the backend falls back to the CPU.
        if (backend != GPUBackend::AUTO) {
            EXPECT_TRUE(gpu.usedBackend() == backend || gpu.usedBackend() == GPUBackend::CPU);
        }
    }
}

// Finds the graph size at which the GPU pays for itself, transfers included.
// NETWORKIT_GPU_BENCHMARK_MAX_SCALE sets the largest graph, with 2^scale nodes and 16 * 2^scale
// edges.
TEST(GPUPageRankBenchmark, compareWithPageRank) {
    Aux::Random::setSeed(42, true);
    const char *maxScaleVariable = std::getenv("NETWORKIT_GPU_BENCHMARK_MAX_SCALE");
    const count maxScale = maxScaleVariable == nullptr ? 22 : std::stoul(maxScaleVariable);

    std::cout << "scale,nodes,edges,cpu_iterations,gpu_iterations,backend,cpu_ms,gpu_ms,cpu/gpu,"
                 "max_error\n";

    constexpr count repetitions = 3;
    for (count scale = 10; scale <= maxScale; scale += 2) {
        const auto G = RmatGenerator(scale, 16, 0.57, 0.19, 0.19, 0.05).generate();

        const auto medianMilliseconds = [&](auto &&algorithm) -> double {
            std::vector<double> times;
            for (count repetition = 0; repetition < repetitions; ++repetition) {
                Aux::Timer timer;
                timer.start();
                algorithm();
                timer.stop();
                times.push_back(static_cast<double>(timer.elapsedMicroseconds()) / 1000.);
            }
            std::sort(times.begin(), times.end());
            return times[times.size() / 2];
        };

        PageRank cpu(G);
        GPUPageRank gpu(G);
        const double cpuMilliseconds = medianMilliseconds([&]() -> void { cpu.run(); });
        const double gpuMilliseconds = medianMilliseconds([&]() -> void { gpu.run(); });

        // Relative to the score: tells how well the ranking survives the GPU's precision.
        double maxError = 0.;
        G.forNodes([&](node u) -> void {
            maxError = std::max(maxError, std::abs(gpu.score(u) - cpu.score(u)) / cpu.score(u));
        });

        std::cout << scale << ',' << G.numberOfNodes() << ',' << G.numberOfEdges() << ','
                  << cpu.numberOfIterations() << ',' << gpu.numberOfIterations() << ','
                  << backendName(gpu.usedBackend()) << ',' << cpuMilliseconds << ','
                  << gpuMilliseconds << ',' << cpuMilliseconds / gpuMilliseconds << ',' << maxError
                  << '\n';
    }
}

} /* namespace NetworKit */
