
#include <algorithm>
#include <cmath>
#include <numeric>

#include <gmock/gmock.h>

#include <girgs/Generator.h>

using namespace std;

// FWD for distance function. Declared in main.
double distance(const std::vector<double>& a, const std::vector<double>& b);

class Generator_test: public testing::Test
{
protected:
    int seed = 1337;
};


bool connected(int a, int b, const vector<pair<int, int>> graph) {
    bool a2b = find(graph.begin(), graph.end(), make_pair(a, b)) != graph.end();
    bool b2a = find(graph.begin(), graph.end(), make_pair(b, a)) != graph.end();
    return a2b || b2a;
}


TEST_F(Generator_test, testThresholdModel)
{
    const auto n = 100;
    const auto alpha = numeric_limits<double>::infinity();
    const auto ple = 2.8;

    auto weights = girgs::generateWeights(n, ple, seed);
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){

        auto positions = girgs::generatePositions(n, d, seed+d);
        auto edges = girgs::generateEdges(weights, positions, alpha, 0);

        // check that there is an edge if and only if the condition in the paper holds: dist < c*(w1w2/W)^-d
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){

                const auto dist = distance(positions[i], positions[j]);
                const auto d_term = pow(dist, d);
                const auto w_term = weights[i] * weights[j] / W;

                if(d_term < w_term) {
                    EXPECT_TRUE(connected(i,j, edges)) << "edge should be present";
                } else {
					EXPECT_FALSE(connected(i,j, edges)) << "edge should be absent";
                }
            }
        }
    }
}

TEST_F(Generator_test, testGeneralModel)
{
    const auto n = 600;
    const auto alpha = 2.5;
    const auto ple = 2.5;

    auto weights = girgs::generateWeights(n, ple, seed);
    auto W = accumulate(weights.begin(), weights.end(), 0.0);

    for(auto d=1u; d<5; ++d){
        // check that the number of generated edges is close to the expected value

        // 1) generator
        auto positions = girgs::generatePositions(n, d, seed+d);
        auto edges = girgs::generateEdges(weights, positions, alpha, seed+d);

        // 2) quadratic sanity check
        auto expectedEdges = 0.0;
        for(int j=0; j<n; ++j){
            for(int i=j+1; i<n; ++i){

                const auto dist = distance(positions[i], positions[j]);
                const auto d_term = pow(dist, d);
                const auto w_term = weights[i] * weights[j] / W;

                auto prob = std::min(std::pow(w_term/d_term, alpha), 1.0);
                expectedEdges += 2*prob;
            }
        }

        auto generatedEdges = edges.size()*2;

        auto rigor = 0.98;
        EXPECT_LT(rigor * expectedEdges, generatedEdges) << "edges too much below expected value";
        EXPECT_LT(rigor * generatedEdges, expectedEdges) << "edges too much above expected value";
    }
}


TEST_F(Generator_test, testCompleteGraph)
{
    const auto n = 100;
    const auto alpha = 0.0; // each edge prob will be 100% now
    const auto ple = 2.5;

    auto weights = girgs::generateWeights(n, ple, seed);

    for(auto d=1u; d<5; ++d) {

        auto positions = girgs::generatePositions(n, d, seed+d);
        auto edges = girgs::generateEdges(weights, positions, alpha, seed+d);

		// check for the correct number of edges
		EXPECT_EQ(edges.size(), (n*(n - 1)) / 2) << "expect a complete graph withour self loops";

        // check that each node is connected to all other nodes
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                EXPECT_TRUE(connected(i,j,edges));
            }
        }
    }
}



// samples all edges by threshold model: dist(i,j) < c*(wiwj/W)^(1/d)
double edgesInQuadraticSampling(const std::vector<double>& w, const vector<vector<double>>& pos, double c) {
    auto n = w.size();
    auto d = pos.front().size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto edges = 0.0;
    for(int i=0; i<n; ++i)
        for(int j=i+1; j<n; ++j)
            if(distance(pos[i], pos[j]) < c*std::pow(w[i] * w[j] / W, 1.0/d))
                edges += 2; // both endpoints get an edge
    return edges;
}


TEST_F(Generator_test, testThresholdEstimation)
{
    auto n = 300;
    auto ple = 2.5;
    auto alpha = numeric_limits<double>::infinity();
    auto weightSeed = seed;
    auto positionSeed = seed;

    auto desired_avg = 10;
    auto runs = 20;

    auto weights = girgs::generateWeights(n, ple, weightSeed);

    // do the tests for all dimensions < 5
    for(auto d = 1; d<5; ++d) {

        // estimate scaling for current dimension
        auto scaled_weights = weights;
        auto scaling = girgs::scaleWeights(scaled_weights, desired_avg, d, alpha);
        auto estimated_c = pow(scaling, 1.0/d);

        // observed avg with estimated c (over multiple runs with different positions)
        auto observed_avg = 0.0;
        for(int i = 0; i<runs; ++i) {

            // try GIRGS generator and quadratic sampling
            auto positions = girgs::generatePositions(n, d, positionSeed+i);
            auto edges = girgs::generateEdges(scaled_weights, positions, alpha, 0);

            auto avg1 = 2.0 * edges.size() / n;
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


TEST_F(Generator_test, testEstimation)
{
    auto all_n = {100, 150, 500};
    auto all_alpha = {0.7, 3.0, numeric_limits<double>::infinity()};
    auto all_desired_avg = {10, 20, 50, 100};
    auto all_dimensions = {1, 2, 3};
    auto runs = 5;

    auto ple = 2.5;
    auto weightSeed = seed;
    auto positionSeed = seed;

    for(int n : all_n){
        for(double alpha : all_alpha){
            for(double desired_avg : all_desired_avg){
                if (desired_avg * 3 > n) continue;
                for(int d : all_dimensions){

                    // generate weights
                    auto weights = girgs::generateWeights(n, ple, weightSeed);

                    // estimate scaling for current dimension
                    girgs::scaleWeights(weights, desired_avg, d, alpha);

                    auto observed_avg = 0.0;
                    for(int i = 0; i<runs; ++i) {

                        // try GIRGS generator
                        auto positions = girgs::generatePositions(n, d, positionSeed+i);
                        auto edges = girgs::generateEdges(weights, positions, alpha, n+i);

                        auto avg = 2.0 * edges.size() / n;
                        observed_avg += avg;
                    }
                    observed_avg /= runs;

                    // test the goodness of the estimation for weight scaling
                    EXPECT_LT(abs(desired_avg - observed_avg)/desired_avg, 0.05) << "estimated constant does not produce desired average degree";
                }
            }
        }
    }
}


TEST_F(Generator_test, testWeightSampling)
{
    auto n = 10000;
    auto ple = 2.1;
    int runs = 10;

    for(int i=0; i<runs; ++i){

        auto weights = girgs::generateWeights(n, ple, seed+i);
        for(auto each : weights) {
            EXPECT_GE(each, 1.0);
            EXPECT_LT(each, n);
        }
        auto max_weight = *max_element(weights.begin(), weights.end());
        EXPECT_GT(max_weight * max_weight, n) << "max weight should be large";
    }
}


TEST_F(Generator_test, testReproducible)
{
    auto n = 1000;
    auto ple = 2.4;
    auto weight_seed    = 1337;
    auto position_seed  = 42;

    auto alphas = { 1.5, std::numeric_limits<double>::infinity() };
    auto dimensions = { 1, 2 };

    for (auto alpha : alphas) {
        for (auto d : dimensions) {

            auto edges1 = girgs::generateEdges(
                    girgs::generateWeights(n, ple, weight_seed),
                    girgs::generatePositions(n, d, position_seed),
                    alpha, weight_seed+position_seed);
            auto edges2 = girgs::generateEdges(
                    girgs::generateWeights(n, ple, weight_seed),
                    girgs::generatePositions(n, d, position_seed),
                    alpha, weight_seed+position_seed);

            // same edges
            EXPECT_EQ(edges1, edges2);
        }
    }
}
