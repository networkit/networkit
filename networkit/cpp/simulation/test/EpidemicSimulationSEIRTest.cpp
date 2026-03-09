/*  EpidemicSimulationSEIRTest.cpp
 *
 *  Created on: 22.02.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gtest/gtest.h>
#include <networkit/graph/Graph.hpp>
#include <networkit/simulation/EpidemicSimulationSEIR.hpp>

#include <array>

namespace NetworKit {

class EpidemicSimulationSEIRTest : public ::testing::Test {
protected:
    static constexpr count states{4};
    static constexpr size_t rowWidth{4}; // {zero, t, s, population}

    using SEIRArray = std::array<count, states>; // {S,E,I,R}

    static count popAt(const std::vector<std::vector<count>> &stats, count maxTimestamps, count t,
                       count s) {
        EXPECT_LT(t, maxTimestamps);
        const count idx = t * states + s; // assumes State::S == 0
        return stats.at(idx).at(3);
    }

    static SEIRArray readSEIR(const std::vector<std::vector<count>> &stats, count maxTimestamps,
                              count t) {
        return {popAt(stats, maxTimestamps, t, 0), popAt(stats, maxTimestamps, t, 1),
                popAt(stats, maxTimestamps, t, 2), popAt(stats, maxTimestamps, t, 3)};
    }

    static void expectTotal(const SEIRArray &SEIR, count n) {
        EXPECT_EQ(SEIR[0] + SEIR[1] + SEIR[2] + SEIR[3], n);
    }

    static void expectStatsShape(const std::vector<std::vector<count>> &stats,
                                 count maxTimestamps) {
        EXPECT_EQ(stats.size(), maxTimestamps * states);
        for (const auto &row : stats) {
            EXPECT_EQ(row.size(), rowWidth);
        }
    }
};

TEST_F(EpidemicSimulationSEIRTest, testConstructor) {
    Graph graph(1);
    constexpr count maxTimestamps{0};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{0};
    constexpr count infectiousTime{0};
    constexpr node startingNode{0};
    EXPECT_NO_THROW(EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb,
                                                     exposureTime, infectiousTime, startingNode));
}

TEST_F(EpidemicSimulationSEIRTest, testNoMaxTimestampsNoStats) {
    Graph graph(2);
    graph.addEdge(0, 1);
    constexpr count maxTimestamps{0};
    constexpr double transmissionProb{1.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{2};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();
    EXPECT_TRUE(stats.empty());
}

TEST_F(EpidemicSimulationSEIRTest, testZeroTransmissionProbabalityNoSpread) {
    Graph graph(10);
    graph.addEdge(0, 1);
    constexpr count maxTimestamps{10};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{2};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    for (count time = 0; time < maxTimestamps; ++time) {
        const auto SEIR = readSEIR(stats, maxTimestamps, time);

        EXPECT_EQ(SEIR[1], 0);
        EXPECT_EQ(SEIR[0], 9);
        if (time >= infectiousTime)
            EXPECT_EQ(SEIR[3], 1);
        else
            EXPECT_EQ(SEIR[2], 1);
        expectTotal(SEIR, 10);
    }
}

TEST_F(EpidemicSimulationSEIRTest, testStatsRowShape) {
    Graph graph(3);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    constexpr count maxTimestamps{4};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{1};
    constexpr count infectiousTime{1};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();
    expectStatsShape(stats, maxTimestamps);
}

TEST_F(EpidemicSimulationSEIRTest, testMaxTimestampsOneRecordsExactlyOneStep) {
    Graph graph(2);
    graph.addEdge(0, 1);
    constexpr count maxTimestamps{1};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{5};
    constexpr count infectiousTime{5};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    for (count time = 0; time < maxTimestamps; ++time) {
        const auto SEIR = readSEIR(stats, maxTimestamps, time);

        EXPECT_EQ(SEIR[0], 1);
        EXPECT_EQ(SEIR[1], 0);
        EXPECT_EQ(SEIR[2], 1);
        EXPECT_EQ(SEIR[3], 0);
        expectTotal(SEIR, 2);
    }
}

TEST_F(EpidemicSimulationSEIRTest, testZeroTransmissionInfectiousTimeZeroImmediateRemoval) {
    Graph graph(10);
    graph.addEdge(0, 1);
    constexpr count maxTimestamps{5};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{0};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    for (count time = 0; time < maxTimestamps; ++time) {
        const auto SEIR = readSEIR(stats, maxTimestamps, time);

        EXPECT_EQ(SEIR[1], 0);
        EXPECT_EQ(SEIR[0], 9);
        EXPECT_EQ(SEIR[2], 0);
        EXPECT_EQ(SEIR[3], 1);
        expectTotal(SEIR, 10);
    }
}

TEST_F(EpidemicSimulationSEIRTest,
       testZeroTransmissionInfectiousTimeGreaterThanMaxTimestampsNeverRemoved) {
    Graph graph(10);
    graph.addEdge(0, 1);
    constexpr count maxTimestamps{5};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{100};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    for (count time = 0; time < maxTimestamps; ++time) {
        const auto SEIR = readSEIR(stats, maxTimestamps, time);

        EXPECT_EQ(SEIR[1], 0);
        EXPECT_EQ(SEIR[0], 9);
        EXPECT_EQ(SEIR[2], 1);
        EXPECT_EQ(SEIR[3], 0);
        expectTotal(SEIR, 10);
    }
}

TEST_F(EpidemicSimulationSEIRTest, testFullTransmissionStarExposureAndRemovalTimeline) {
    Graph starGraph(6);
    for (node v = 1; v < 6; ++v) {
        starGraph.addEdge(0, v);
    }

    constexpr count maxTimestamps{5};
    constexpr double transmissionProb{1.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{2};
    constexpr node startingNode{0};

    EpidemicSimulationSEIR simulator(starGraph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    for (count time = 0; time < maxTimestamps; ++time) {
        const auto SEIR = readSEIR(stats, maxTimestamps, time);

        if (time == 0 || time == 1) {
            EXPECT_EQ(SEIR[0], 0);
            EXPECT_EQ(SEIR[1], 5);
            EXPECT_EQ(SEIR[2], 1);
            EXPECT_EQ(SEIR[3], 0);
        } else if (time == 2 || time == 3) {
            EXPECT_EQ(SEIR[0], 0);
            EXPECT_EQ(SEIR[1], 0);
            EXPECT_EQ(SEIR[2], 5);
            EXPECT_EQ(SEIR[3], 1);
        } else if (time == 4) {
            EXPECT_EQ(SEIR[0], 0);
            EXPECT_EQ(SEIR[1], 0);
            EXPECT_EQ(SEIR[2], 0);
            EXPECT_EQ(SEIR[3], 6);
        }

        expectTotal(SEIR, 6);
    }
}

TEST_F(EpidemicSimulationSEIRTest,
       testFullTransmissionStarWithLongExposureTimeNoSecondaryInfections) {
    Graph starGraph(6);
    for (node v = 1; v < 6; ++v) {
        starGraph.addEdge(0, v);
    }

    constexpr count maxTimestamps{4};
    constexpr double transmissionProb{1.0};
    constexpr count exposureTime{100};
    constexpr count infectiousTime{2};
    constexpr node startingNode{0};

    EpidemicSimulationSEIR simulator(starGraph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    for (count time = 0; time < maxTimestamps; ++time) {
        const auto SEIR = readSEIR(stats, maxTimestamps, time);

        if (time < infectiousTime) {
            EXPECT_EQ(SEIR[0], 0);
            EXPECT_EQ(SEIR[1], 5);
            EXPECT_EQ(SEIR[2], 1);
            EXPECT_EQ(SEIR[3], 0);
        } else {
            EXPECT_EQ(SEIR[0], 0);
            EXPECT_EQ(SEIR[1], 5);
            EXPECT_EQ(SEIR[2], 0);
            EXPECT_EQ(SEIR[3], 1);
        }

        expectTotal(SEIR, 6);
    }
}

TEST_F(EpidemicSimulationSEIRTest, testStartNoneChoosesRandomNode) {
    Graph graph(4);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);

    constexpr count maxTimestamps{1};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{2};
    constexpr node startingNode{none};

    EpidemicSimulationSEIR simulator(graph, maxTimestamps, transmissionProb, exposureTime,
                                     infectiousTime, startingNode);
    simulator.run();
    auto stats = simulator.getData();

    expectStatsShape(stats, maxTimestamps);

    const count chosenStart = stats.front().at(0);
    EXPECT_NE(chosenStart, none);
    EXPECT_LT(chosenStart, graph.upperNodeIdBound());

    for (const auto &row : stats) {
        EXPECT_EQ(row.at(0), chosenStart);
    }

    const auto SEIR = readSEIR(stats, maxTimestamps, 0);
    EXPECT_EQ(SEIR[0], 3);
    EXPECT_EQ(SEIR[1], 0);
    EXPECT_EQ(SEIR[2], 1);
    EXPECT_EQ(SEIR[3], 0);
    expectTotal(SEIR, 4);
}
} // namespace NetworKit
