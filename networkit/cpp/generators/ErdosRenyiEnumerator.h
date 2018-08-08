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
public:
	ErdosRenyiEnumerator(node n, double prob, bool directed) :
		n{n},
		prob{prob},
		inv_log_cp{1.0 / std::log(1.0 - prob)},
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


private:
	const node n;
	const double prob;
	const double inv_log_cp;
	const bool directed;

// Directed is easier as we simply iterated over the complete
// adjacency matrix
	template <typename Handle>
	void enumerateDirected(Handle handle, unsigned tid, const node node_begin, const node node_end) const {
		Aux::SignalHandler handler;

	// random source
		auto& prng = Aux::Random::getURNG(); // this is thread local
		std::uniform_real_distribution<double> distr{
			std::nextafter(0.0, 1.0), std::nextafter(1.0, 0.0)};

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
 		std::uniform_real_distribution<double> distr{
			std::nextafter(0.0, 1.0), std::nextafter(1.0, 0.0)};

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
	count skip_distance(double random_prob) const {
		assert(0 < random_prob && random_prob < 1);
		return 1 + static_cast<count>(floor(log(random_prob) * inv_log_cp));
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