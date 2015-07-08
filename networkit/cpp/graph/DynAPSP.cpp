/*
 * DynAPSP.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#include "APSP.h"
#include "DynAPSP.h"
// for testing
#include "Dijkstra.h"
// end for testing
#include "../auxiliary/Log.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/NumericTools.h"
#include <queue>


namespace NetworKit {

DynAPSP::DynAPSP(const Graph& G) : APSP(G) {

}

void DynAPSP::dynamic_sssp(node root, node startbfs) {
	distances[0][0] = 42;
}

void DynAPSP::update(const GraphEvent& event) {
	if (event.type!=GraphEvent::EDGE_ADDITION && event.type!=GraphEvent::EDGE_WEIGHT_UPDATE)
		throw std::runtime_error("Graph update not allowed");
	// we assume edge insertion
	// node u = event.u;
	// node v = event.v;
	//
	// std::vector<std::vector<edgeweight> > L;
	// L = distances;
	// G.forNodes([&](node x){
	// 	if (L[x][v] - L[x][u] > 1) {
	// 		L[x][v] = distances[x][u] + 1;
	// 		INFO("x, v, L[x, v]: ", x," ", v," ", L[x][v]);
	// 		dynamic_sssp(x, v);
	// 	}
	// });

	// we assume edge deletion
	node u = event.u;
	node v = event.v;

	std::vector<std::vector<edgeweight> > L;
	L = distances;

	// update l.[v] block
	Graph Gtrans = G.transpose();
	// use dijkstra for now. We might be able to use dyn_sssp_1
	Dijkstra dijk(Gtrans, v);
	dijk.run();
	std::vector<edgeweight> distancesToV = dijk.getDistances();
	G.forNodes([&](node x){
		L[x][v] = distancesToV[x];
	// end of update l.[v] block
		if (L[x][v] - L[x][u] > 1) {
			INFO("x, v, L[x, v]: ", x," ", v," ", L[x][v]);
			dynamic_sssp(x, v);
		}
	});

	// mod = false;
	// // priority queue with distance-node pairs
	// Aux::PrioQueue<edgeweight, node> Q(G.upperNodeIdBound());
	// // queue with all visited nodes
	// std::queue<node> visited;
	// // if u has a new shortest path going through v, it updates the distance of u
	// // and inserts u in the priority queue (or updates its priority, if already in Q)
	// auto updateQueue = [&](node u, node v, edgeweight w) {
	// 	if (distances[u] >= distances[v]+w) {
	// 		distances[u] = distances[v]+w;
	// 		if (color[u] == WHITE) {
	// 			Q.insert(distances[u], u);
	// 			color[u] = BLACK;
	// 		}	else {
	// 			Q.decreaseKey(distances[u], u);
	// 		}
	// 	}
	// };
	//
	// for (GraphEvent edge : batch) {
	// 	if (edge.type!=GraphEvent::EDGE_ADDITION && edge.type!=GraphEvent::EDGE_WEIGHT_UPDATE)
	// 		throw std::runtime_error("Graph update not allowed");
	// 	//TODO: discuss with Christian whether you can substitute weight_update with with_increase/weight_decrease
	// 	// otherwise, it is not possbile to check wether the change in the weight is positive or negative
	// 	updateQueue(edge.u, edge.v, edge.w);
	// 	updateQueue(edge.v, edge.u, edge.w);
	// }
	//
	// while(Q.size() != 0) {
	// 	mod = true;
	// 	node current = Q.extractMin().second;
	// 	visited.push(current);
	// 	if (storePreds) {
	// 		previous[current].clear();
	// 	}
	// 	npaths[current] = 0;
	// 	G.forInNeighborsOf(current, [&](node current, node z, edgeweight w){
	// 		//z is a predecessor of current node
	// 		if (Aux::NumericTools::equal(distances[current], distances[z]+w, 0.000001)) {
	// 			if (storePreds) {
	// 				previous[current].push_back(z);
	// 			}
	// 			npaths[current] += npaths[z];
	// 		}
	// 		//check whether curent node is a predecessor of z
	// 		else {
	// 			updateQueue(z, current, w);
	// 		}
	// 	});
	// }
	// reset colors
	// while(!visited.empty()) {
	// 	node w = visited.front();
	// 	visited.pop();
	// 	color[w] = WHITE;
	// }
}

} /* namespace NetworKit */
