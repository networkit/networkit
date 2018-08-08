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


class ErdosRenyiEnumerator {
	//! this type is used only internally for fixed-point arithmetic
	using integral_t = unsigned long long;

public:
	ErdosRenyiEnumerator(node n, double prob, bool directed) :
		n{n},
		inv_log2_cp{1.0 / std::log2(1.0 - prob)},
		directed{directed}
	{
		assert(n > 0);
		assert(0 < prob);
		assert(prob < 1);
	}

	template<typename Handle>
	void forEdges(Handle handle) {
		if (directed) {
			enumerateDirected(handle, 0, 0, n);
		} else {
			enumerateUndirected(handle, 0, 0, n);
		}
	}

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

				enumerateDirected(handle, tid, first_node, last_node);
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
					enumerateUndirected(handle, tid, first_node, last_node);
			}
		}
	}

	/**
	 * Similarly to @a forEdgesParallel but computed only on one thread.
	 * If the callback accepts three arguments tid is always 0.
	 */
	template<typename Handle>
	void forEdges(Handle handle) {
		if (directed) {
			enumerateDirected(handle, 0, 0, n);
		} else {
			enumerateUndirected(handle, 0, 0, n);
		}
	}


private:
	const node n;
	const double inv_log2_cp;
	const bool directed;

// Directed is easier as we simply iterated over the complete
// adjacency matrix
	template <typename Handle>
	void enumerateDirected(Handle handle, unsigned tid, const node node_begin, const node node_end) const {
		Aux::SignalHandler handler;

	// random source
		auto& prng = Aux::Random::getURNG(); // this is thread local
		std::uniform_int_distribution<integral_t> distr{1, std::numeric_limits<integral_t>::max()};

		count curr = node_begin;
		node next = -1;
		while (curr < node_end) {
			handler.assureRunning();
			next += skip_distance(distr(prng));

			curr += next / n;
			next = next % n;

			if (curr < node_end) {
				callHandle(handle, tid, curr, next);
			}

 		}
	}

// In the undirected case we only traverse the lower triangle (excluding the
// diagonal) of the adjacency matrix
	template <typename Handle>
	void enumerateUndirected(Handle handle, unsigned tid, const node node_begin, const node node_end) const {
		Aux::SignalHandler handler;

		// random source
		auto& prng = Aux::Random::getURNG(); // this is thread local
 		std::uniform_int_distribution<integral_t> distr{1, std::numeric_limits<integral_t>::max()};

		count curr = std::max<node>(node_begin, 1u);
		node next = -1;

		while (curr < node_end) {
			handler.assureRunning();
			// compute new step length
			next += skip_distance(distr(prng));

			// check if at end of row
			while (next >= curr) {
				// adapt to next row
				next = next - curr;
				curr++;
			}

			// insert edge
			if (curr < node_end) {
				callHandle(handle, tid, curr, next);
			}
		}
	}

	// Skip Magic due to Batagelj and Brandes
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