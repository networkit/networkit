/*
 * TriangleCounter.cpp
 *
 *  Created on: 29.08.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include "TriangleCounter.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"
#include <omp.h>

namespace NetworKit {

std::vector<int> TriangleCounter::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	auto isOutEdge = [&](node u, node v) {
		return graph.degree(u) > graph.degree(v) || (graph.degree(u) == graph.degree(v) && u < v);
	};

	std::vector<std::vector<std::pair<node, edgeid> > > inEdges(graph.upperNodeIdBound());

	graph.balancedParallelForNodes([&](node u) {
		// bucket sort
		std::vector<count> nodePos(1);

		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid){
			if (isOutEdge(v, u)) {
				count deg = graph.degree(v);
				if (nodePos.size() < deg + 1) {
					nodePos.resize(deg + 1);
				}

				++nodePos[deg];
			}
		});
		
		// bucket sort
		// exclusive prefix sum
		index tmp = nodePos[0];
		index sum = tmp;
		nodePos[0] = 0;
		for (index i = 1; i < nodePos.size(); ++i) {
			tmp = nodePos[i];
			nodePos[i] = sum;
			sum += tmp;
		}

		inEdges[u].resize(sum);

		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid){
			if (isOutEdge(v, u)) {
				inEdges[u][nodePos[graph.degree(v)]++] = std::make_pair(v, eid);
			}
		});
	});

	//Edge attribute: triangle count
	std::vector<int> triangleCount(graph.upperEdgeIdBound(), 0);

	std::vector<std::vector<edgeid> > nodeMarker;

	#pragma omp parallel
	{
		#pragma omp single
		{
			nodeMarker.resize(omp_get_num_threads());
		}

		nodeMarker[omp_get_thread_num()].resize(graph.upperNodeIdBound(), none);
	}

	graph.balancedParallelForNodes([&](node u) {
		auto tid = omp_get_thread_num();
		count degU = graph.degree(u);
		std::vector<std::pair<node, edgeid> > outEdges;
		outEdges.reserve(degU - inEdges[u].size());

		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			if (isOutEdge(u, v)) {
				outEdges.emplace_back(v, eid);
				nodeMarker[tid][v] = 0;
			}
		});

		//For all neighbors: check for already marked neighbors.
		for (auto uv : outEdges) {
			for (auto vw : inEdges[uv.first]) {
				if (graph.degree(vw.first) > degU) break;

				if (nodeMarker[tid][vw.first] != none) {

					++nodeMarker[tid][uv.first];
					++nodeMarker[tid][vw.first];
					#pragma omp atomic update
					++triangleCount[vw.second];
				}
			}
		}

		for (auto uv : outEdges) {
			if (nodeMarker[tid][uv.first] > 0) {
				#pragma omp atomic update
				triangleCount[uv.second] += nodeMarker[tid][uv.first];
			}
			nodeMarker[tid][uv.first] = none;
		}
	});

	return triangleCount;
}

} /* namespace NetworKit */
