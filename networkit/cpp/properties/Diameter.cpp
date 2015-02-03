/*
 * Diameter.cpp
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#include <numeric>

#include "Diameter.h"
#include "Eccentricity.h"
#include "../graph/BFS.h"
#include "../graph/Dijkstra.h"
#include "../properties/ConnectedComponents.h"
#include "../graph/BFS.h"

namespace NetworKit {


edgeweight Diameter::exactDiameter(const Graph& G) {
	using namespace std;

	edgeweight diameter = 0.0;

	if (! G.isWeighted()) {
		std::tie(diameter, std::ignore) = estimatedDiameterRange(G, 0);
	} else {
		 G.forNodes([&](node v) {
		 	Dijkstra dijkstra(G, v);
		 	dijkstra.run();
		 	auto distances = dijkstra.getDistances();
		 	for (auto distance : distances) {
		 		if (diameter < distance) {
		 			diameter = distance;
		 		}
		 	}
//			DEBUG("ecc(", v, "): ", *std::max_element(distances.begin(), distances.end()), " of ", distances);
		 });
	}

	if (diameter == std::numeric_limits<edgeweight>::max()) {
		throw std::runtime_error("Graph not connected - diameter is infinite");
	}
	return diameter;
}





std::pair<edgeweight, edgeweight> Diameter::estimatedDiameterRange(const NetworKit::Graph &G, double error, std::pair< NetworKit::node, NetworKit::node > *proof) {
	if (G.isDirected()) {
		throw std::runtime_error("Error, the diameter of directed graphs cannot be computed yet.");
	}
	/*
	 * This is an implementation of the iFub-algorithm by
	 * Pilu Crescenzi, Roberto Grossi, Michel Habib, Leonardo Lanzi, Andrea Marino:
	 * On computing the diameter of real-world undirected graphs,
	 * Theoretical Computer Science, Volume 514, 25 November 2013, Pages 84-95, ISSN 0304-3975,
	 * http://dx.doi.org/10.1016/j.tcs.2012.09.018.
	 * (http://www.sciencedirect.com/science/article/pii/S0304397512008687)
	 */
	// start node selection. Here the very simple variant with the node of max. degree
	// for each component
	std::vector<node> startNodes;
	{
		ConnectedComponents comp(G);
		comp.run();
		count numberOfComponents = comp.numberOfComponents();
		startNodes.resize(numberOfComponents, 0);
		std::vector<count> maxDeg(numberOfComponents, 0);

		// for each component, find the node with the maximum degreee and add it as start node
		G.forNodes([&](node v) {
			count d = G.degree(v);
			count c = comp.componentOfNode(v);
			if (d > maxDeg[c]) {
				startNodes[c] = v;
				maxDeg[c] = d;
			}
		});
	}

	// simple BFS from all start nodes, store all levels (needed later)
	std::vector<std::vector<node> > level(1);
	G.BFSfrom(startNodes, [&](node v, count dist) {
		if (level.size() < dist+1) level.emplace_back();
		level[dist].push_back(v);
	});
	count eccStart = level.size()-1;

	// set initial lower and upper bound
	count lb = eccStart, ub = 2*eccStart;

	// calculate the ecc for the nodes in each level starting with the highest level
	// until either the bound is tight or all nodes have been considered
	for (count i = eccStart; ub > (lb + error * lb) && i > 0; --i) {
		// we need to wait until here as otherwise after a break in the inner loop the wrong ub might
		// be returned if it was decreased too early
		ub = 2*i;
		INFO("New upper bound ", ub);
		// second possibility, if we did not meet the break condition we might now meet it
		if (ub <= (lb + error * lb)) {
			break;
		}

		bool aborted = false;

		#pragma omp parallel for schedule(guided)
		for (index j = 0; j < level[i].size(); ++j) {
			#pragma omp flush (aborted)
			if (!aborted) {
				node v = level[i][j];
				count ecc = 0;
				node w = none;
				std::tie(w, ecc) = Eccentricity::getValue(G, v);
				// check if the lower bound has been improved and if yes, if the upper bound has been met
				if (ecc > lb) {
					#pragma omp critical
					if (ecc > lb) {
						lb = ecc;

						INFO("Found new lower bound ", lb);

						if (proof != NULL) {
							proof->first = v;
							proof->second = w;
						}

						if (ub <= (lb + error * lb)) {
							aborted = true;
						}
					}
				}
			}
		}
	}

	return {lb, std::max(ub, lb)};
}


edgeweight Diameter::estimatedVertexDiameter(const Graph& G, count samples) {

	edgeweight infDist = std::numeric_limits<edgeweight>::max();

	// TODO: consider weights

	auto estimateFrom = [&](node v) -> count {
		BFS bfs(G, v);
		bfs.run();
		auto distances = bfs.getDistances();

		// get two largest path lengths
		edgeweight maxD = 0;
		edgeweight maxD2 = 0; // second largest distance
		for (auto d : distances) {
			if ((d != infDist) && (d >= maxD)) {
				maxD2 = maxD;
				maxD = d;
			}
		}

		edgeweight dMax = maxD + maxD2;
		count vd = (count) dMax + 1; 	// count the nodes, not the edges
		return vd;
	};

	edgeweight vdMax = 0;
	#pragma omp parallel for
	for (count i = 0; i < samples; ++i) {
		node u = G.randomNode();
		edgeweight vd = estimateFrom(u);
		DEBUG("sampled vertex diameter from node ", u, ": ", vd);
		#pragma omp critical
		{
			if (vd > vdMax) {
				vdMax = vd;
			}
		}
	}

	return vdMax;

}


edgeweight Diameter::estimatedVertexDiameterPedantic(const Graph& G) {

	edgeweight infDist = std::numeric_limits<edgeweight>::max();

	auto estimateFrom = [&](node v) -> count {
		BFS bfs(G, v);
		bfs.run();
		auto distances = bfs.getDistances();

		// get two largest path lengths
		count maxD = 0;
		count maxD2 = 0; // second largest distance
		for (auto d : distances) {
			if ((d != infDist) && (d >= maxD)) {
				maxD2 = maxD;
				maxD = d;
			}
		}

		edgeweight dMax = maxD + maxD2;
		count vd = (count) dMax + 1; 	// count the nodes, not the edges
		return vd;
	};

	ConnectedComponents cc(G);
	DEBUG("finding connected components");
	cc.run();
	INFO("Graph is weighted? ", G.isWeighted());
	if (!G.isWeighted()) {
		INFO("Number of components ", cc.numberOfComponents());
		if (cc.numberOfComponents() > 1) {
			DEBUG("estimating for each component in parallel");
			std::vector<std::set<node> > components;
			for (auto component : cc.getPartition().getSubsets()) {
				components.push_back(component);
			}
			DEBUG("gathered components");
			std::vector<count> vds;
			#pragma omp parallel for
			for (index i = 0; i < components.size(); ++i) {
				count vd = estimateFrom(*components[i].begin()); // take any node from the component and perform bfs from there
				DEBUG("checking component ", i);
				#pragma omp critical
				vds.push_back(vd);
			}

			count vdMax = *std::max_element(vds.begin(), vds.end());
			return vdMax;

		} else {
			return estimateFrom(G.randomNode());
		}
	} else {
		INFO("Number of components ", cc.numberOfComponents());
		DEBUG("Estimating size of the largest component");
		std::map<count, count> sizes = cc.getComponentSizes();
		count largest_comp_size = 0;
		for(auto it = sizes.cbegin(); it != sizes.cend(); ++it) {
			DEBUG(it->second);
			if (it->second > largest_comp_size) {
				largest_comp_size = it->second;
			}
		}
		INFO("Largest component size: ", largest_comp_size);
		return largest_comp_size;
	}

}


} /* namespace NetworKit */
