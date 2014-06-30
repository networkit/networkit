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
		G.forNodes([&](node v) {
			BFS bfs(G, v);
			bfs.run();
			auto distances = bfs.getDistances();
			for (auto distance : distances) {
				if (diameter < distance) {
					diameter = distance;
				}
			}

//			DEBUG("ecc(", v, "): ", *std::max_element(distances.begin(), distances.end()), " of ", distances);
		});
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





std::pair<edgeweight, edgeweight> Diameter::estimatedDiameterRange(const Graph& G, double error) {
	/* BFS that calls f with the visited edges and returns the node with largest distance from u. */
	/* Note: the function Graph::breadthFirstEdgesFrom that should
	 do the same has not been implemented!
		-- TODO: Then why not implement it directly there?
	 */
	auto bfs_edges = [&] (const Graph& G, node u, std::function<void(node, node)> f) -> node {
		std::queue<node> q;
		std::vector<bool> visited(G.upperNodeIdBound(), false);
		q.push(u);
		visited[u] = true;

		node x = u;
		while (!q.empty()) {
			x = q.front(); q.pop();
			G.forNeighborsOf(x, [&] (node y) {
						if (!visited[y]) {
							f(x, y);
							visited[y] = true;
							q.push(y);
						}
					});
		}
		return x;
	};

	/* Diameter estimate: lowerBounds <= diam(G) <= upperBound. */
	edgeweight lowerBound = 0.0;
	edgeweight upperBound = std::numeric_limits <edgeweight>::max();
	const count z = G.upperNodeIdBound();

	/* Nodes sorted decreasingly by degree. */
	std::vector<node> high_deg(z);
	std::iota(begin(high_deg), end(high_deg), 0);
	std::sort(begin(high_deg), end(high_deg), [&] (node u, node v) {
		// if (G.hasNode(u) && G.hasNode(v)) { // TODO: what if nodes are not present
		return G.degree(u) > G.degree(v);
	});

	/* Random node. */
	// TODO: use random node routine in Graph
	static const std::default_random_engine random;
	auto random_node = std::bind(std::uniform_int_distribution<node>(0, z - 1), random);

	/* While not converged: update estimate. */
	count niter = 0;
	while ((upperBound - lowerBound) >= error * lowerBound && niter < G.numberOfNodes()) {
		edgeweight ecc;

		/* ecc(u) <= diam(G) */
		node u = random_node();
//		DEBUG("u: ", u);
		std::tie(std::ignore, ecc) = Eccentricity::getValue(G, u);
		lowerBound = std::max(lowerBound, ecc);

		/* diam(G) <= diam(BFS_Tree(v)) */
		node v = high_deg[niter];
		Graph bfs_tree(z);
		node w = bfs_edges(G, v, [&] (node a, node b) {
			bfs_tree.addEdge(a, b);
		});
//		DEBUG("v: ", v, ", w: ", w);

		bfs_tree.forEdges([&](node u, node v) {
			assert(G.hasEdge(u, v));
		});


		/* diam(T) = ecc_T(w) by problem 4. */
		std::tie(std::ignore, ecc) = Eccentricity::getValue(bfs_tree, w);
		upperBound = std::min(upperBound, ecc);

		niter++;
	}

	if ((lowerBound == std::numeric_limits<edgeweight>::max()) || (upperBound == std::numeric_limits<edgeweight>::max())) {
		throw std::runtime_error("Graph not connected - diameter is infinite");
	}

	return {lowerBound, upperBound};
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

}


} /* namespace NetworKit */
