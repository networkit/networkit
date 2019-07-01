#include <algorithm>
#include <cmath>
#include <numeric>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include "../Generator.hpp"

class GirgsGeneratorGTest: public testing::Test {
public:
    static double distance(const std::vector<double>& a, const std::vector<double>& b) {
        assert(a.size() == b.size());
        auto result = 0.0;
        for(auto d=0u; d<a.size(); ++d){
            auto dist = std::abs(a[d] - b[d]);
            dist = std::min(dist, 1.0-dist);
            result = std::max(result, dist);
        }
        return result;
    }
};

TEST_F(GirgsGeneratorGTest, testThresholdModel)
{
    const auto n = 100;
    constexpr auto alpha = std::numeric_limits<double>::infinity();
    const auto ple = 2.8;

    Aux::Random::setSeed(1332, true);

    auto weights = NetworKit::girgs::generateWeights(n, ple);
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){

        auto positions = NetworKit::girgs::generatePositions(n, d);
        auto graph = NetworKit::girgs::generateEdges(weights, positions, alpha);

        // check that there is an edge if and only if the condition in the paper holds: dist < c*(w1w2/W)^-d
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){

                const auto dist = GirgsGeneratorGTest::distance(positions[i], positions[j]);
                const auto d_term = pow(dist, d);
                const auto w_term = weights[i] * weights[j] / W;

                if(d_term < w_term) {
                    EXPECT_TRUE(graph.hasEdge(i,j)) << "edge should be present";
                } else {
					EXPECT_FALSE(graph.hasEdge(i,j)) << "edge should be absent";
                }
            }
        }
    }
}

TEST_F(GirgsGeneratorGTest, testGeneralModel)
{
    const auto n = 600;
    const auto alpha = 2.5;
    const auto ple = 2.5;

    Aux::Random::setSeed(1333, true);

    auto weights = NetworKit::girgs::generateWeights(n, ple);
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){
        // check that the number of generated edges is close to the expected value

        // 1) generator
        auto positions = NetworKit::girgs::generatePositions(n, d);
        auto graph = NetworKit::girgs::generateEdges(weights, positions, alpha);

        // 2) quadratic sanity check
        auto expectedEdges = 0.0;
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){

                const auto dist = GirgsGeneratorGTest::distance(positions[i], positions[j]);
                const auto d_term = pow(dist, d);
                const auto w_term = weights[i] * weights[j] / W;

                auto prob = std::min(std::pow(w_term/d_term, alpha), 1.0);
                expectedEdges += prob;
            }
        }

        auto generatedEdges = graph.numberOfEdges();

        auto rigor = 0.98;
        EXPECT_LT(rigor * expectedEdges, generatedEdges) << "edges too much below expected value";
        EXPECT_LT(rigor * generatedEdges, expectedEdges) << "edges too much above expected value";
    }
}


TEST_F(GirgsGeneratorGTest, testCompleteGraph)
{
    const auto n = 100;
    const auto alpha = 0.0; // each edge prob will be 100% now
    const auto ple = 2.5;

    Aux::Random::setSeed(1334, true);
    auto weights = NetworKit::girgs::generateWeights(n, ple);

    for(auto d=1u; d<5; ++d) {

        auto positions = NetworKit::girgs::generatePositions(n, d);
        auto graph = NetworKit::girgs::generateEdges(weights, positions, alpha);

		// check for the correct number of edges
		EXPECT_EQ(graph.numberOfEdges(), (n*(n - 1)) / 2) << "expect a complete graph withour self loops";

        // check that each node is connected to all other nodes
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                EXPECT_TRUE(graph.hasEdge(i,j));
            }
        }
    }
}

// samples all edges by threshold model: dist(i,j) < c*(wiwj/W)^(1/d)
static double edgesInQuadraticSampling(const std::vector<double>& w, const std::vector<std::vector<double>>& pos, double c) {
    auto n = w.size();
    auto d = pos.front().size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto edges = 0.0;
    for(unsigned i=0; i<n; ++i)
        for(unsigned j=i+1; j<n; ++j)
            if(GirgsGeneratorGTest::distance(pos[i], pos[j]) < c*std::pow(w[i] * w[j] / W, 1.0/d))
                edges += 2; // both endpoints get an edge
    return edges;
}


TEST_F(GirgsGeneratorGTest, testThresholdEstimation)
{
    const auto n = 300;
    const auto ple = 2.5;
    constexpr auto alpha = std::numeric_limits<double>::infinity();

    auto desired_avg = 10;
    auto runs = 20;

    Aux::Random::setSeed(1336, true);
    auto weights = NetworKit::girgs::generateWeights(n, ple);

    // do the tests for all dimensions < 5
    for(auto d = 1; d<5; ++d) {

        // estimate scaling for current dimension
        auto scaled_weights = weights;
        auto scaling = NetworKit::girgs::scaleWeights(scaled_weights, desired_avg, d, alpha);
        auto estimated_c = pow(scaling, 1.0/d);

        // observed avg with estimated c (over multiple runs with different positions)
        auto observed_avg = 0.0;
        for(int i = 0; i<runs; ++i) {

            // try GIRGS generator and quadratic sampling
            auto positions = NetworKit::girgs::generatePositions(n, d);
            auto graph = NetworKit::girgs::generateEdges(scaled_weights, positions, alpha);

            auto avg1 = 2.0 * graph.numberOfEdges() / n;
            auto avg2 = edgesInQuadraticSampling(weights, positions, estimated_c) / n;

            // generator must yield same results as quadratic sampling
            EXPECT_EQ(avg1, avg2) << "sampling with scaled weights produced different results than quadratic samping with constant factor";
            observed_avg += avg1;
        }
        observed_avg /= runs;

        // test the goodness of the estimation for weight scaling
        EXPECT_LT(abs(desired_avg - observed_avg) / desired_avg, 0.05) << "estimated constant does not produce desired average degree";
    }
}


TEST_F(GirgsGeneratorGTest, testEstimation)
{
    auto all_n = {100, 150, 250};
    auto all_alpha = {0.7, 3.0, std::numeric_limits<double>::infinity()};
    auto all_desired_avg = {10, 20, 100};
    auto all_dimensions = {1, 2, 3};
    auto runs = 10;

    auto ple = 2.5;

    Aux::Random::setSeed(1335, true);
    for(int n : all_n) {
        for(double alpha : all_alpha) {
            for(double desired_avg : all_desired_avg) {
                if (desired_avg * 3 > n) continue;
                for(int d : all_dimensions){

                    // generate weights
                    auto weights = NetworKit::girgs::generateWeights(n, ple);

                    // estimate scaling for current dimension
                    NetworKit::girgs::scaleWeights(weights, desired_avg, d, alpha);

                    auto observed_avg = 0.0;
                    for(int i = 0; i<runs; ++i) {
                        // try GIRGS generator
                        auto positions = NetworKit::girgs::generatePositions(n, d);
                        auto graph = NetworKit::girgs::generateEdges(weights, positions, alpha);

                        observed_avg += 2*graph.numberOfEdges();
                    }
                    observed_avg /= runs;
                    observed_avg /= n;

                    // test the goodness of the estimation for weight scaling
					EXPECT_LT(abs(desired_avg - observed_avg) / desired_avg, 0.05)
						<< " d=" << d << " n=" << n << " alpha=" << alpha 
						<< " desired_avg=" << desired_avg << " observed_avg=" << observed_avg;
                }
            }
        }
    }
}


TEST_F(GirgsGeneratorGTest, testWeightSampling)
{
    auto n = 10000;
    auto ple = 2.1;
    int runs = 10;

    Aux::Random::setSeed(1337, true);
    for(int i=0; i<runs; ++i){
        auto weights = NetworKit::girgs::generateWeights(n, ple);
        for(auto each : weights) {
            EXPECT_GE(each, 1.0);
            EXPECT_LT(each, n);
        }
        auto max_weight = *max_element(weights.begin(), weights.end());
        EXPECT_GT(max_weight * max_weight, n) << "max weight should be large";
    }
}


TEST_F(GirgsGeneratorGTest, testReproducible)
{
    auto n = 1000;
    auto ple = 2.4;

    auto alphas = { 1.5, std::numeric_limits<double>::infinity() };
    auto dimensions = { 1, 2 };

    for (auto alpha : alphas) {
        for (auto d : dimensions) {
            auto get_graph = [&] {
                Aux::Random::setSeed(1337, true);

                auto weights = NetworKit::girgs::generateWeights(n, ple);
                auto positions = NetworKit::girgs::generatePositions(n, d);

                return NetworKit::girgs::generateEdges(weights, positions, alpha);
            };

            auto graph1 = get_graph();
            auto graph2 = get_graph();

            // there seems to be no equality check for NetworKit's graphs
            EXPECT_EQ(graph1.numberOfEdges(), graph2.numberOfEdges());
            graph1.forEdges([&graph2] (NetworKit::node a, NetworKit::node b) {
                EXPECT_TRUE(graph2.hasEdge(a,b));
            });
        }
    }
}
