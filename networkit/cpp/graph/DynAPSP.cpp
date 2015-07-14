/*
 * DynAPSP.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
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
	INFO("Entering Dynamic SSSP, root = ", root, ", startbfs = ", startbfs);
	// priority queue with distance-node pairs
	Aux::PrioQueue<edgeweight, node> Q(G.upperNodeIdBound());
	std::vector<bool> enqueued(G.upperNodeIdBound(), false);
	// queue with all visited nodes
	// if u has a new shortest path going through v, it updates the distance of u
	// and inserts u in the priority queue (or updates its priority, if already in Q)
	auto updateQueue = [&](node u, edgeweight priority) {
		if (!enqueued[u]) {
			Q.insert(priority, u);
			enqueued[u] = true;
		}	else {
			Q.decreaseKey(priority, u);
		}
	};

	updateQueue(startbfs, std::min(L[root][startbfs], distances[root][startbfs]));
	while(Q.size() != 0) {
		//INFO("Entering while loop");
		node v = Q.extractMin().second;
		INFO("The current node v is ", v);
		enqueued[v] = false;
		INFO("label[",root,"][",v,"] = ", L[root][v]," distances[",root,"][",v,"] = ", distances[root][v]);
		if (L[root][v] < distances[root][v]) {
			INFO("label[",root,"][",v,"] < distances[",root,"][",v,"]");
			INFO("distances[",root,"][",v,"] was ", distances[root][v], " and now becomes ",L[root][v]);
			distances[root][v] = L[root][v];
			if(!G.isDirected()) {
				distances[v][root] = L[root][v];
			}
			G.forNeighborsOf(v, [&](node w, edgeweight weight){
				INFO("neighbor w: ", w);
				if (distances[root][v] + weight < L[root][w]) {
					L[root][w] = distances[root][v] + weight;
					updateQueue(w, std::min(L[root][w], distances[root][w]));
				}
			});
		}
		if (L[root][v] > distances[root][v]) {
			//INFO("label[",root,"][",v,"] > distances[",root,"][",v,"]");
			edgeweight Dold = distances[root][v];
			distances[root][v] = infDist;
			edgeweight con = infDist;
			G.forInNeighborsOf(v, [&](node z, edgeweight weight) {
				if (distances[root][z] != infDist && distances[root][z] + weight < con) {
					con = distances[root][z] + weight;
				}
			});
			L[root][v] = con;
			updateQueue(v, L[root][v]);
			G.forNeighborsOf(v, [&](node w, edgeweight weight){
				if(Dold + weight == L[root][w]) {
					edgeweight con = infDist;
					G.forInNeighborsOf(w, [&](node z, edgeweight weight) {
						if (distances[root][z] != infDist && distances[root][z] + weight < con) {
							con = distances[root][z] + weight;
						}
					});
					L[root][w] = con;
					updateQueue(w, std::min(L[root][w], distances[root][w]));
				}
			});
		}
	}
}

void DynAPSP::update(const GraphEvent& event) {
	INFO("Entering update");
	if (event.type!=GraphEvent::EDGE_ADDITION && event.type!=GraphEvent::EDGE_REMOVAL && event.type!=GraphEvent::EDGE_WEIGHT_UPDATE)
		throw std::runtime_error("Graph update not allowed");

	node u = event.u;
	node v = event.v;
	L = distances;

	if (event.type==GraphEvent::EDGE_ADDITION) {
		G.forNodes([&](node x){
			if (L[x][v] - L[x][u] > 1) {
				L[x][v] = distances[x][u] + 1;
				INFO("x, v, D[x][v], L[x, v]: ", x, " ", v," ", distances[x][v], " ", L[x][v]);
				dynamic_sssp(x, v);
			}
		});
	}
	if (event.type==GraphEvent::EDGE_REMOVAL) {
		// update l.[v] block
		// use dijkstra for now. We might be able to use dyn_sssp_1
		std::vector<edgeweight> distancesToV(G.upperNodeIdBound());
		if (G.isDirected()) {
			Graph Gtrans = G.transpose();
			Dijkstra dijk(Gtrans, v);
			dijk.run();
			distancesToV = dijk.getDistances();
		} else {
			Dijkstra dijk(G, v);
			dijk.run();
			distancesToV = dijk.getDistances();
		}
		
		INFO("v, distancesToV:", v, " ", distancesToV[1]," ", distancesToV[2]," ", distancesToV[3]," ", distancesToV[4]," ", distancesToV[5]," ", distancesToV[6]);
		G.forNodes([&](node x){
			L[x][v] = distancesToV[x];
		// end of update l.[v] block
			if (L[x][v] - L[x][u] > 1) {
				INFO("x, v, D[x][v], L[x, v]: ", x, " ", v," ", distances[x][v], " ", L[x][v]);
				dynamic_sssp(x, v);
			}
		});
	}
}

} /* namespace NetworKit */
