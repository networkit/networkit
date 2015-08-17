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

void DynAPSP::dynamic_sssp(node root, std::vector<std::pair<node, edgeweight> > batch) {
	INFO("Entering Dynamic SSSP, root = ", root);
	INFO("Labels: ", L[root][0], " ", L[root][1], " ", L[root][2], " ", L[root][3], " ", L[root][4]);
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

	for (std::pair<node, edgeweight> update: batch) {
		INFO("Enqueueing node ", update.first, " with priority ", update.second);
		updateQueue(update.first, update.second);
	}

	while(Q.size() != 0) {
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
				if (distances[root][v] + weight < L[root][w]) {
					L[root][w] = distances[root][v] + weight;
					INFO("neighbor ", w, " of node ", v, " inserted into queue.");
					updateQueue(w, std::min(L[root][w], distances[root][w]));
				}
			});
		}
		if (L[root][v] > distances[root][v]) {
			INFO("label[",root,"][",v,"] > distances[",root,"][",v,"]");
			edgeweight Dold = distances[root][v];
			distances[root][v] = infDist;
			edgeweight con = infDist;
			G.forInNeighborsOf(v, [&](node z, edgeweight weight) {
				if (distances[root][z] != infDist && distances[root][z] + weight < con) {
					con = distances[root][z] + weight;
				}
			});
			L[root][v] = con;
			INFO("reinserting ", v, " with priority ", L[root][v]);
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
					INFO("inserting neighbor ", w," of ", v, " with priority ", L[root][v]);
					updateQueue(w, std::min(L[root][w], distances[root][w]));
				}
			});
		}
	}
}

/**
 * update takes a graph event and updates the all-pair shortest path distances. L contains
 * the label value for each node. A node's label value is necessary for the sssp
 * recomputation (see Tuned SWSF in Bauer & Wagner paper titled Batch Dynamic Sinlge-Source
 * Shortest-Path Algorithms: An Experimental Study). update initializes the label values and
 * identifies all nodes for which an sssp must be run. The helper method dynamic_sssp does the
 * actual updating of the distances.
 */
void DynAPSP::update(const std::vector<GraphEvent>& batch) {
	INFO("Entering update");
	// initialize affected sources
	std::map<node, std::vector<std::pair<node, edgeweight> > > affected;
	std::map<node, std::vector<std::pair<node, edgeweight> > >::iterator affectedIterator;

	L = distances;
	
	// identify affected sources
	for (GraphEvent event: batch) {
		if (event.type!=GraphEvent::EDGE_ADDITION && event.type!=GraphEvent::EDGE_REMOVAL &&
			event.type!=GraphEvent::EDGE_WEIGHT_UPDATE && event.type!=GraphEvent::EDGE_WEIGHT_INCREMENT)
			throw std::runtime_error("Graph update not allowed");

		node u = event.u;
		node v = event.v;
		edgeweight ew = G.weight(u, v);

		if (event.type==GraphEvent::EDGE_ADDITION || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && ew < 0)
			|| (event.type==GraphEvent::EDGE_WEIGHT_UPDATE && ew < G.weight(u,v)) ) {
			G.forNodes([&](node x){
				if (L[x][v] > ew + distances[x][u]) {
					INFO("edge insertion. Updating L[x][v]. x, v, old, new : ", x, " ", v, " ", L[x][v], " ", distances[x][u] + ew);
					L[x][v] = distances[x][u] + ew;
					INFO("x, v, D[x][v], L[x, v]: ", x, " ", v," ", distances[x][v], " ", L[x][v]);
					affected[x].push_back(std::pair<node, edgeweight>(v, std::min(L[x][v], distances[x][v])));
				}
			});
		}
		if (event.type==GraphEvent::EDGE_REMOVAL || (event.type==GraphEvent::EDGE_WEIGHT_INCREMENT && ew > 0)
			|| (event.type==GraphEvent::EDGE_WEIGHT_UPDATE && ew > G.weight(u,v)) ) {

			// this part below identifies all source nodes for which v must be set to con(v)

			G.forNodes([&](node x){
				edgeweight con = infDist;
				if (x != v) {
					G.forInNeighborsOf(v, [&](node y, edgeweight weightyv) {
						if (distances[x][y] != infDist && distances[x][y] + weightyv < con) {
							con = distances[x][y] + weightyv;
						}
					});
					INFO("edge removal. Setting L[x][v] to con. x, v, L[x][v], con: ", x, " ", v, " ", L[x][v], " ", con);
					L[x][v] = con;
				}

				//INFO("x, v, L[x][u], L[x][v], L[x][v] - L[x][u], ew: ", x, " ", v," ", L[x][u], " ", L[x][v], " ", L[x][v] - L[x][u], " ", ew);
				if (L[x][v] > ew + distances[x][u]) {
					INFO("x, v, D[x][v], L[x, v]: ", x, " ", v," ", distances[x][v], " ", L[x][v]);
					affected[x].push_back(std::pair<node, edgeweight>(v, std::min(L[x][v], distances[x][v])));
				}
			});
		}
	}

	// process affected sources
	for (affectedIterator = affected.begin(); affectedIterator != affected.end(); affectedIterator++) {
		dynamic_sssp(affectedIterator->first, affectedIterator->second);
	}
}

} /* namespace NetworKit */
