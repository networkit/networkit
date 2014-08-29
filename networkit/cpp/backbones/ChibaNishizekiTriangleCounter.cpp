/*
 * ChibaNishizekiTriangleCounter.cpp
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#include "ChibaNishizekiTriangleCounter.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

std::vector<int> ChibaNishizekiTriangleCounter::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	std::vector<std::vector<std::pair<node, edgeid> > > edges(graph.upperNodeIdBound());
	
	graph.parallelForNodes([&](node u) {
		edges[u].reserve(graph.degree(u));
		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			edges[u].emplace_back(v, eid);
		});
	});

	//Node attribute: marker
	std::vector<edgeid> nodeMarker(graph.upperNodeIdBound(), none);

	//Edge attribute: triangle count
	std::vector<int> triangleCount(graph.upperEdgeIdBound(), 0);

	graph.forNodes([&](node u) {
		//Mark all neighbors
		for (auto uv : edges[u]) {
			nodeMarker[uv.first] = uv.second;
		}

		//For all neighbors: check for already marked neighbors.
		for (auto uv : edges[u]) {
			bool edgeDeleted = false;
			for (auto vw = edges[uv.first].begin(); vw != edges[uv.first].end(); ++vw) {
				if (edgeDeleted) {
					(*(vw-1)) = *vw;
				}
				if (nodeMarker[vw->first] != none) {

					edgeid eid_uw = nodeMarker[vw->first];

					++triangleCount[uv.second];
					++triangleCount[eid_uw];
					++triangleCount[vw->second];
				} else if (vw->first == u) {
					edgeDeleted = true;
				}
			}
			
			assert(edgeDeleted);
			
			edges[uv.first].pop_back();

			nodeMarker[uv.first] = none;
		}
	});

	return triangleCount;
}

} /* namespace NetworKit */
