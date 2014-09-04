/*
 * TopDegreeAttributizer.cpp
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#include "TopDegreeAttributizer.h"
#include "LocalSimilarityAttributizer.h" 	//For AttributizedEdge
#include <limits>

#include "../auxiliary/Log.h"

namespace NetworKit {

TopDegreeAttributizer::TopDegreeAttributizer() {}

std::vector<count> TopDegreeAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {

	std::vector<count> sparsificationAttribute (graph.upperEdgeIdBound(),  std::numeric_limits<count>::max());
	graph.forNodes([&](node i) {
		//Get a list of neighbors, in descending order by node degree
		std::vector<AttributizedEdge<count>> neighbors;
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			neighbors.push_back(AttributizedEdge<count> (i, j, eid, graph.degree(j)));
		});
		std::sort(neighbors.begin(), neighbors.end(),
			[](const AttributizedEdge<count>& n1, const AttributizedEdge<count>& n2) {
					return n1 < n2; } );

		count rank = 1;
		count previousValue = std::numeric_limits<count>::max();
		for (std::vector<AttributizedEdge<count>>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
			edgeid eid = it-> eid;

			sparsificationAttribute[eid] = std::min(rank, sparsificationAttribute[eid]);

			if (previousValue != it->value) {
				previousValue = it->value;
				rank++;
			}
		}

	});

	return sparsificationAttribute;
}

} /* namespace NetworKit */
