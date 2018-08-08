#include "ErdosRenyiEnumeratorGTest.h"
#include "../ErdosRenyiEnumerator.h"

#include "../../auxiliary/Random.h"

namespace NetworKit {

static inline size_t edgeindex(node u, node v, node n) {
	return u*n + v;
}

TEST_P(ErdosRenyiEnumeratorGTest, TestSequential) {
	const bool directed = std::get<0>(GetParam());
	const node n        = std::get<1>(GetParam());
	const double prob   = std::get<2>(GetParam());

	Aux::Random::setSeed(static_cast<size_t>((123456 + directed * 23425) * prob * n), false);

	std::vector<unsigned> edge_counts(n * n, 0);
	std::vector<bool> edge_exists(n * n, 0);

	const edgeid expected_edges = static_cast<edgeid>(
		(directed ? n*n : (n*n - n) / 2)  * prob);

	// after this many rounds each edge should be hit at least once
	const int rounds = 10 * (1/prob) * std::log(1/prob);

	for(int i=0; i<rounds; i++) {
		ErdosRenyiEnumerator ere(n, prob, directed);

		// initialize data
		edge_exists.assign(n*n, false);
		count num_edges = 0;

		ere.forEdges([&] (node u, node v) {
			// ensure edge is legal (for undirected edges the
			// first node is larger than the second)
			ASSERT_LE(u, n);
			ASSERT_LE(v, directed ? n : u);

			// check that edge is unique
			ASSERT_FALSE(edge_exists[edgeindex(u,v,n)]);
			edge_exists[edgeindex(u,v,n)] = true;

			// collect statistics for probabilisitic checks
			edge_counts[edgeindex(u,v,n)]++;
			num_edges++;
		});

		// count that we produced a number of nodes, near the expected number
		// we could apply much sharper bounds here
		ASSERT_GE(num_edges, expected_edges/2);
		ASSERT_LE(num_edges, 2*expected_edges);
	}

	for(node u = 1; u != n; u++) {
		const node upper = directed ? n : u;
		for(node v = 0; v != upper; v++) {
			// Again, much sharper bounds are possible if we do not
			// rely on the implicit union bound here.
			auto count = edge_counts[edgeindex(u,v,n)];
			EXPECT_GE(count, 1) << "edge(" << u << ", " << v << "), rounds" << rounds;
			ASSERT_LE(count, log(rounds) * rounds * prob) << "edge(" << u << ", " << v << "), rounds" << rounds;
		}
	}
}


TEST_P(ErdosRenyiEnumeratorGTest, TestParallel) {
	const bool directed = std::get<0>(GetParam());
	const node n        = std::get<1>(GetParam());
	const double prob   = std::get<2>(GetParam());

	Aux::Random::setSeed(static_cast<size_t>((1234567 + directed * 23425) * prob * n), true);

	std::vector<unsigned> edge_counts(n * n, 0);
	std::vector<unsigned> edge_exists(n * n, 0); // std::vector<bool> is not thread safe

	const edgeid expected_edges = static_cast<edgeid>(
		(directed ? n*n : (n*n - n) / 2)  * prob);

	// after this many rounds each edge should be hit at least once
	const int rounds = 10 * (1/prob) * std::log(1/prob);

	for(int i=0; i<rounds; i++) {
		ErdosRenyiEnumerator ere(n, prob, directed);

		// initialize data
		edge_exists.assign(n*n, false);
		constexpr size_t cacheline_scale = 17; // avoid false sharing
		std::vector<count> num_edges_thread(omp_get_max_threads() * cacheline_scale);

		ere.forEdgesParallel([&] (int tid, node u, node v) {
			// ensure edge is legal (for undirected edges the
			// first node is larger than the second)
			ASSERT_LE(u, n);
			ASSERT_LE(v, directed ? n : u);

			// check that edge is unique
			ASSERT_FALSE(edge_exists[edgeindex(u,v,n)]);
			edge_exists[edgeindex(u,v,n)] = true;

			// collect statistics for probabilisitic checks
			edge_counts[edgeindex(u,v,n)]++;
			num_edges_thread[cacheline_scale*tid]++;
		});

		size_t num_edges = std::accumulate(num_edges_thread.cbegin(), num_edges_thread.cend(), 0u);
		for(auto count : num_edges_thread)
			ASSERT_LT(count, num_edges);

		// count that we produced a number of nodes, near the expected number
		// we could apply much sharper bounds here
		ASSERT_GE(num_edges, expected_edges/2);
		ASSERT_LE(num_edges, 2*expected_edges);

	}

	for(node u = 1; u != n; u++) {
		const node upper = directed ? n : u;
		for(node v = 0; v != upper; v++) {
			// Again, much sharper bounds are possible if we do not
			// rely on the implicit union bound here.
			auto count = edge_counts[edgeindex(u,v,n)];
			EXPECT_GE(count, 1) << "edge(" << u << ", " << v << "), rounds" << rounds;
			ASSERT_LE(count, log(rounds) * rounds * prob) << "edge(" << u << ", " << v << "), rounds" << rounds;
		}
	}
}



INSTANTIATE_TEST_CASE_P(ErdosRenyiEnumeratorGTest, ErdosRenyiEnumeratorGTest,
						::testing::Values(
						 std::make_tuple(false, 100, 0.1),    std::make_tuple(true, 100, 0.1),
						 std::make_tuple(false, 200, 0.005),  std::make_tuple(true, 200, 0.005)
						));

} // ! namespace NetworKit