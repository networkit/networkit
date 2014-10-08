/*
 *
 */

#include "ChibaNishizekiQuadrangleCounter.h"

std::vector< int > NetworKit::ChibaNishizekiQuadrangleCounter::getAttribute(const NetworKit::Graph &graph, const std::vector< int > &attribute) {
	std::vector<std::vector<std::pair<node, edgeid> > > edges(graph.upperNodeIdBound());

	// copy edges with edge ids
	graph.parallelForNodes([&](node u) {
		edges[u].reserve(graph.degree(u));
		graph.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			edges[u].emplace_back(v, eid);
		});
	});

	//Node attribute: marker
	std::vector<count> nodeMarker(graph.upperNodeIdBound(), 0);

	//Edge attribute: triangle count
	std::vector<int> quandrangleCount(graph.upperEdgeIdBound(), 0);

	// bucket sort
	count n = graph.numberOfNodes();
	std::vector<node> sortedNodes(n);
	{
		std::vector<index> nodePos(n + 1, 0);

		graph.forNodes([&](node u) {
			++nodePos[n - graph.degree(u)];
		});

		// exclusive prefix sum
		index tmp = nodePos[0];
		index sum = tmp;
		nodePos[0] = 0;

		for (index i = 1; i < nodePos.size(); ++i) {
			tmp = nodePos[i];
			nodePos[i] = sum;
			sum += tmp;
		}

		graph.forNodes([&](node u) {
			sortedNodes[nodePos[n - graph.degree(u)]++] = u;
		});
	}

	for (node u : sortedNodes) {
		for (auto uv : edges[u]) {
			for (auto vw = edges[uv.first].begin(); vw != edges[uv.first].end(); ++vw) {
				// delete the edge to u as we do not need to consider it again.
				// the opposite edge doesn't need to be deleted as we will never again consider
				// outgoing edges of u as u cannot be reached anymore after the uv loop.
				if (vw->first == u) {
					// move last element to current position in order to avoid changing too much
					*vw = edges[uv.first].back();
					edges[uv.first].pop_back();
					if (vw == edges[uv.first].end()) // break if we were at the last element already
						break;
				}

				++nodeMarker[vw->first];
			}
		}

		for (auto uv : edges[u]) {
			for (auto vw : edges[uv.first]) {
				if (nodeMarker[vw.first] > 1) {
					quandrangleCount[uv.second] += nodeMarker[vw.first] - 1;
					quandrangleCount[vw.second] += nodeMarker[vw.first] - 1;
				}
			}
		}

		for (auto uv : edges[u]) {
			for (auto vw : edges[uv.first]) {
				nodeMarker[vw.first] = 0;
			}
		}
	}

	return quandrangleCount;

}
