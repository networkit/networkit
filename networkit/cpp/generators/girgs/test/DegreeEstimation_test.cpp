
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include <girgs/Generator.h>


using namespace std;

// FWD for distance function. Declared in main.
double distance(const std::vector<double>& a, const std::vector<double>& b);

// multiple functions follow to compute the expected number of edges
double tryOften(const std::vector<double>& w, double c, int dim, double a) {

    auto n = w.size();

    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto gen = std::mt19937(random_device()());
    std::uniform_real_distribution<> dist; // [0..1)

    int runs = 20;
    double avg = 0.0;
    for(int k=0; k<runs; ++k){

        auto sum = 0.0;
        // sample different positions for each run
        auto pos = vector<vector<double>>(n, vector<double>(dim, 0.0));
        for(int i=0; i<n; ++i)
            for (int d=0; d<dim; ++d)
                pos[i][d] = dist(gen);

        for(int i=0; i<n; ++i)
            for(int j=i+1; j<n; ++j) {
                auto w_term = w[i] * w[j] / W;
                auto d_term = pow(distance(pos[i], pos[j]), dim);
                auto edgeProb = min(c * pow(w_term/d_term, a), 1.0);
                sum += 2*edgeProb;
            }
        avg += sum;
    }
    avg /= runs;
    return avg;
}


double shortEdgesTrivial(const std::vector<double>& w, double c, int d, double a) {
    auto n = w.size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);

    auto short_edges = 0.0;
    for(int i=0; i<n; ++i)
        for(int j=i+1; j<n; ++j) {
            auto short_edge_prob = (1<<d) * pow(c, 1/a) * w[i]*w[j]/W;
            short_edges += 2* min(short_edge_prob, 1.0); // MINIMUM IS MISSING IN THE THESIS
        }
    return short_edges;
}

double longEdgesTrivial(const std::vector<double>& w, double c, int d, double a) {
    auto n = w.size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);

    // long edges
    auto long_edges = 0.0;
    for(int i=0; i<n; ++i)
        for(int j=i+1; j<n; ++j) {
            auto w_term = w[i] * w[j] / W;
            auto crazy_w = pow(c, 1/a/d) * pow(w_term, 1.0/d);

            if(crazy_w >= 0.5) // INTERVAL BOUNDS ARE WRONG IF w > 1/2
                continue;

            // since crazy_w is the turning point for the minimum we must get Pr = 1 for dist = crazy_W
            auto control_one = c * pow(w_term / pow(crazy_w, d), a);
            EXPECT_LT(0.9999, control_one);
            EXPECT_LT(control_one, 1.0001);

            auto long_edge_prob = c*pow(w_term, a) * d*(1<<d) / (d - a*d) * ( pow(0.5, d-a*d) - pow(crazy_w, d-a*d));
            EXPECT_LE(long_edge_prob, 1.0); // prob for a long edge cannot be bigger than 1

            long_edges += 2*long_edge_prob;
        }

    return long_edges;
}


// takes w by value so not to sort the original weights
vector<double> getRichClub(std::vector<double> w, double c, int d, double a) {

    auto n = w.size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);

    vector<double> rich_club;
    sort(w.begin(), w.end(), std::greater<double>());
    auto w_n = w.front();
    for(int i=0; i<n; ++i){
        auto crazy_w = pow(c, 1/a/d) * pow(w[i] * w_n / W, 1.0/d);
        if(crazy_w > 0.5)
            rich_club.push_back(w[i]);
        else
            break;
    }
    return rich_club;
}


double shortEdgesImproved(const std::vector<double>& w, double c, int d, double a) {

    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto rich_club = getRichClub(w,c,d,a);

    // compute short edges improved
    auto sum_sq_w   = 0.0; // sum_{v\in V} (w_v^2/W)
    for(auto each : w)
        sum_sq_w   += each*each/W;
    auto shortEdgesWithError = (1<<d) * pow(c, 1/a) * (W-sum_sq_w);
    auto error = 0.0;
    for(int i = 0; i<rich_club.size(); ++i)
        for(int j = 0; j<rich_club.size(); ++j) {
            if(i==j) continue;
            auto w1 = rich_club[i];
            auto w2 = rich_club[j];
            auto e = max( (1<<d)*pow(c, 1/a)*(w1*w2/W) -1.0, 0.0);
            error += e;
            if(e <= 0) break;
        }

    return shortEdgesWithError - error;
}


double longEdgesImproved(const std::vector<double>& w, double c, int d, double a) {

    auto n = w.size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);

    // long edges with error quadratic
    auto long_edges = 0.0;
    for(int i=0; i<n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            auto w_term = w[i] * w[j] / W;
            auto crazy_w = pow(c, 1 / a / d) * pow(w_term, 1.0 / d);
            auto long_edge_prob = c * pow(w_term, a) * d * (1 << d) / (d - a * d) * (pow(0.5, d - a * d) - pow(crazy_w, d - a * d));
            long_edges += 2 * long_edge_prob;
        }
    }
    auto longEdgesWithError = long_edges;

    auto rich_club = getRichClub(w,c,d,a);

    auto error = 0.0;
    for(int i = 0; i<rich_club.size(); ++i)
        for(int j = 0; j<rich_club.size(); ++j) {
            if(i==j) continue;
            auto w_term = rich_club[i] * rich_club[j] / W;
            auto crazy_w = pow(c, 1 / a / d) * pow(w_term, 1.0 / d);
            if(crazy_w <= 0.5)
                break;

            auto e = c * pow(w_term, a) * d * (1 << d) / (d - a * d) * (pow(0.5, d - a * d) - pow(crazy_w, d - a * d));
            error += e;
        }

    return longEdgesWithError - error;
}



double finalForm(const std::vector<double>& w, double c, int d, double a){

    auto n = w.size();
    auto W = std::accumulate(w.begin(), w.end(), 0.0);
    auto sum_sq_w   = 0.0; // sum_{v\in V} (w_v^2/W)
    auto sum_w_a    = 0.0; // sum_{v\in V} (w_v  /W)^\alpha
    auto sum_sq_w_a = 0.0; // sum_{v\in V} (w_v^2/W)^\alpha
    for(auto each : w){
        sum_sq_w   += each*each/W;
        sum_w_a    += pow(each/W, a);
        sum_sq_w_a += pow(each*each/W, a);
    }

    //   sum_{u\in V} sum_{v\in V} (wu*wv/W)^\alpha
    // = sum_{u\in V} sum_{v\in V} wu^\alpha * (wv/W)^\alpha
    // = sum_{u\in V} wu^\alpha sum_{v\in V} (wv/W)^\alpha
    auto sum_wwW_a = 0.0;
    for(auto each : w)
        sum_wwW_a += pow(each, a)*sum_w_a;

    auto factor1 = (W-sum_sq_w) * (1+1/(a-1)) * (1<<d);
    auto factor2 = pow(2, a*d) / (a-1) * (sum_wwW_a - sum_sq_w_a);

    auto long_and_short_with_error = pow(c, 1/a) * factor1 - c * factor2;


    // get error for long and short edges
    auto rich_club = getRichClub(w,c,d,a);
    auto short_error = 0.0;
    auto long_error = 0.0;
    for(int i = 0; i<rich_club.size(); ++i) {
        for (int j = 0; j < rich_club.size(); ++j) {
            if (i == j) continue;

            auto w_term = rich_club[i] * rich_club[j] / W;
            auto crazy_w = pow(c, 1 / a / d) * pow(w_term, 1.0 / d);
            if (crazy_w <= 0.5)
                break;

            short_error += (1<<d)*pow(c, 1/a)*w_term -1.0;
            long_error += c * pow(w_term, a) * d * (1 << d) / (d - a * d) * (pow(0.5, d - a * d) - pow(crazy_w, d - a * d));
        }
    }

    return long_and_short_with_error - short_error - long_error;
}



TEST(DegreeEstimation_test, testEstimationFormula)
{
    const auto seed = 42;
    auto n = 300;
    auto a = 4.5;
    auto d = 2;
    auto ple = 2.5;

    auto c = 0.5;

    auto epsilon = 0.00001;

    auto weights = girgs::generateWeights(n, ple, seed);

    auto experimental_number_of_edges = tryOften(weights, c, d, a);

    auto short_edges = shortEdgesTrivial(weights,c,d,a);
    auto long_edges = longEdgesTrivial(weights,c,d,a);
    auto edges = short_edges + long_edges;

    auto test_short = shortEdgesImproved(weights,c,d,a);
    auto test_long = longEdgesImproved(weights,c,d,a);
    auto test_edges = finalForm(weights, c,d,a);

    EXPECT_LT(abs(short_edges - test_short), epsilon);
    EXPECT_LT(abs(long_edges  - test_long ), epsilon);
    EXPECT_LT(abs(edges       - test_edges), epsilon);

    EXPECT_LT(0.99*edges, experimental_number_of_edges);
    EXPECT_GT(1.01*edges, experimental_number_of_edges);
}
