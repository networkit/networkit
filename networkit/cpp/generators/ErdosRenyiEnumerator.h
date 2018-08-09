/*
 * ErdosRenyiGenerator.h
 *
 *  Created on: 07.08.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef ERDOSRENYIENUMERATOR_H_
#define ERDOSRENYIENUMERATOR_H_

#include <omp.h>

#include <cassert>
#include <cmath>
#include <random>

#include "../Globals.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/SignalHandling.h"

namespace NetworKit {

/**
 * Generates a stream of edges of a G(n,p) graph. The edges are not
 * written to memory, but emitted only via usual @a forEdges semantics
 * to a callback.
 *
 * Use @ref ErdosRenyiGenerator as a wrapper to output a graph.
 *
 * The enumerator implements both floating point and fixed point arithmetic,
 * to compute the edges which can be selected via the template parameter
 * @a UseFixedPoint. It defaults to true, as this algorithm is typically
 * 2 to 3 times faster. In theory, floating point arithmetic allows for
 * longer runs of consecutive edges not selected. As those events virtually
 * never occur, there are no measureable implications of using the faster
 * variant.
 */
template <bool UseFixedPoint = true>
class ErdosRenyiEnumerator {
	//! this type is used only internally for fixed-point arithmetic
	using integral_t = unsigned long long;

public:
	/**
	 * Generates a G(n, p) graph for n > 1 and 0 < p < 1.
	 *
	 * For an @b directed graph, the resulting edge stream is equivalent
	 * to throwing a coin for each node pair (u, v) and accepting it
	 * with probability @a prob; i.e. the stream may include self-loops.
	 * Hence the expected number of edges is n*n*prob.
	 *
	 * For an @b undirected graph, all node pairs (u, v) with 1 < v < u < n
	 * are considered. Hence the expected number of edges is n*(n-1)/2*prob
	 *
	 * @param n			Number of nodes to generate
	 * @param prob      Probability that an edge exists
	 * @param directed  Selects an directed graph
	 */
	ErdosRenyiEnumerator(node n, double prob, bool directed) :
		n{n},
		prob{prob},
		inv_log2_cp{1.0 / std::log2(1.0 - prob)},
		directed{directed}
	{
		assert(n > 0);
		assert(0 < prob);
		assert(prob < 1);
	}

	/**
	 * Generates an Erdos-Renyi-Graph as specified in the constructor on the fly.
	 * The stream is generated in parallel on as many core as an OpenMP parallel
	 * section is allotted to. For each edge the callback @a handle is invoked.
	 * Two signatures are supported for callback:
	 *
	 *   (unsigned tid, node u, node v)
	 *   (node u, node v)
	 *
	 * where tid is the current thread id as returned by omp_get_thread_num() and
	 * u and v are the edge's nodes. In case of an undirected graph u > v.
	 *
	 * It is guaranteed that no two threads emit edges for the same u.
	 *
	 * It can be expected that all threads emit a similar number of edges.
	 */
	template<typename Handle>
	void forEdgesParallel(Handle handle) {
		if (directed) {
			#pragma omp parallel
			{
				const unsigned threads = omp_get_num_threads();
				const unsigned tid = omp_get_thread_num();

				const node chunk_size = (n + threads - 1) / threads;
				const node first_node = std::min<node>(n, tid * chunk_size);
				const node last_node  = std::min<node>(n, (tid+1) * chunk_size);

				enumerate<true>(handle, tid, first_node, last_node);
			}
		} else {
			#pragma omp parallel
			{
				const unsigned threads = omp_get_num_threads();
				const unsigned tid = omp_get_thread_num();

				// cells in adj matrix per thread
				const node chunk_size = (n * (n+1) / 2 + threads - 1) / threads;

				node first_node;
				node last_node = 0;

				for(unsigned i = 0; i < tid; i++) {
					first_node = last_node;
					node nodes_in_chunk = std::ceil(std::sqrt(
						0.5 + first_node * first_node + first_node + 2*chunk_size));
					last_node = std::min<node>(n, first_node + nodes_in_chunk);
				}

				if (first_node < last_node)
					enumerate<false>(handle, tid, first_node, last_node);
			}
		}
	}

	/**
	 * Similarly to @ref forEdgesParallel but computed on one thread only.
	 * If the callback accepts three arguments tid is always 0.
	 */
	template<typename Handle>
	void forEdges(Handle handle) {
		if (directed) {
			enumerate<true>(handle, 0, 0, n);
		} else {
			enumerate<false>(handle, 0, 0, n);
		}
	}

	/**
	 * Returns the expected number of edges to be generated.
	 */
	count expectedNumberOfEdges() const {
		return prob * n * (directed ? n : (n-1) * 0.5);
	}


private:
	const node n; //< number of nodes
	const double prob; //< probability p
	const double inv_log2_cp; //<  1.0 / log2(1 - prob)
	const bool directed; //< true if a directed graph should be generated

// In the undirected case we only traverse the lower triangle (excluding the
// diagonal) of the adjacency matrix
	template <bool Directed, typename Handle>
	void enumerate(Handle handle, unsigned tid, const node node_begin, const node node_end) const {
		Aux::SignalHandler handler;

		// random source
		auto& prng = Aux::Random::getURNG(); // this is thread local
		auto distr = get_distribution<UseFixedPoint>();

		count curr = node_begin;
		if (!Directed && !curr) curr = 1;
		node next = -1;

		node max_skip = 0;

		while (curr < node_end) {
			handler.assureRunning();
			// compute new step length
			auto skip = skip_distance(distr(prng));
			next += skip;
			if (skip > max_skip) max_skip = skip;

			// check if at end of row (assuming an average degree of 1,
			// its typically faster to repeatedly subtract than to
			// divide; it is in any case easier easier ...)
			while (next >= (Directed ? n : curr)) {
				// adapt to next row
				next = next - (Directed ? n : curr);
				curr++;
			}

			// insert edge
			if (curr < node_end) {
				callHandle(handle, tid, curr, next);
			}
		}

		std::cout << max_skip << "\n";
	}

// Optimized version of the computation of the skip distance as
// proposed Batagelj and Brandes. It basically converts a uniform
// variate to a binomial
	count skip_distance(integral_t random_prob) const {
		/*
		 * The original idea is to compute
		 * 1 + floor(log(1 - x) / log(1 - prob)) where x is
		 * a uniform real from (0, 1).
		 *
		 * We now precompute inv_log_cp := 1.0 / log(1 - prob) once
		 * to avoid recomputing the second log and to replace a
		 * division by a multiplication. Hence we compute
		 * 1 + floor(log(1 - x) * inv_log_cp).
		 *
		 * Observe that P[x = k] = P[x = 1-k] and hence we avoid the subtraction:
		 * 1 + floor(log(x) * inv_log_cp).
		 *
		 * Then, typically log() is slower than log2(). On my
		 * machines its a factor of roughly 1.5. Thus we replace
		 * log by log2 in both equations:
		 *  inv_log2_cp := 1.0 / log2(1 - prob)
		 *  1 + floor(log2(x) * inv_log2_cp)
		 *
		 * Now we optimise the generation of the random number.
		 * uniform_real_distribution is roughly 3 times slower than
		 * uniform_int_distribution. Hence let's switch to fix-point arithmetic.
		 * Let X a real drawn uniformly from (0, 2**64), i.e. we model
		 * X = (2**64) * x:
		 *
		 *    1 + floor(log2(x) * inv_log2_cp)
		 *  = 1 + floor(log2(X/2**64) * inv_log2_cp)
		 *  = 1 + floor((log2(X) - 64) * inv_log2_cp).
		 */
		return 1 + static_cast<count>(
			floor((log2(random_prob) - 8*sizeof(integral_t)) * inv_log2_cp)
	    );
	}

	count skip_distance(double random_prob) const {
		return 1 + static_cast<count>(floor((log2(random_prob)) * inv_log2_cp));
	}

// SFINAE to determine and construct the right uniform distribution
	template <bool FixedPoint>
	auto get_distribution() const -> typename std::enable_if<FixedPoint, std::uniform_int_distribution<integral_t>>::type {
		return std::uniform_int_distribution<integral_t>{1, std::numeric_limits<integral_t>::max()};
	}

	template <bool FixedPoint>
	auto get_distribution() const -> typename std::enable_if<!FixedPoint, std::uniform_real_distribution<double>>::type {
		return std::uniform_real_distribution<double>{std::nextafter(0.0, 1.0), std::nextafter(1.0, 0.0)};
	}

// SFINAE to allow using functors with and without thread id as parameter
	template <typename Handle>
	auto callHandle(Handle h, unsigned tid, node u, node v) const -> decltype(h(0u, node{0}, node{0})) {
		return h(tid, u, v);
	}

	template <typename Handle>
	auto callHandle(Handle h, unsigned /*tid*/, node u, node v) const -> decltype(h(node{0}, node{0})) {
		return h(u, v);
	}
};

}

#endif // ERDOSRENYIENUMERATOR_H_