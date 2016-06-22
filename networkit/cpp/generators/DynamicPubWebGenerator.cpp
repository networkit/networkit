/*
 * DynamicPubWebGenerator.cpp
 *
 *  Created on: 15.01.2014
 *      Author: Henning
 */

#include <queue>

#include "DynamicPubWebGenerator.h"

namespace NetworKit {

DynamicPubWebGenerator::DynamicPubWebGenerator(count numNodes,
		count numberOfDenseAreas, float neighborhoodRadius,
		count maxNumberOfNeighbors, bool writeInitialGraphToStream) :
		initGen(numNodes, numberOfDenseAreas, neighborhoodRadius,
				maxNumberOfNeighbors), writeInitialGraphToStream(writeInitialGraphToStream)
{
	G = initGen.generate();
}

std::vector<GraphEvent> DynamicPubWebGenerator::generate(count nSteps) {
	count numToDel = (count) (G.numberOfNodes() * 0.05); // TODO: externalize, possibly randomize
	count numToIns = (count) (G.numberOfNodes() * 0.05); // TODO: externalize, possibly randomize
	std::vector<GraphEvent> eventStream;
	std::vector<std::pair<node, node> > edgesToDelete;
	coordinates.clear();

	if (writeInitialGraphToStream) {
		// write initial graph to stream
		G.forNodes([&](node v){
			eventStream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, v));
		});

		G.forEdges([&](node u, node v, edgeweight ew){
			eventStream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, u, v, ew));
		});	

		eventStream.push_back(GraphEvent(GraphEvent::TIME_STEP));

		writeInitialGraphToStream = false;
	}


	for (index step = 0; step < nSteps; ++step) {

		// delete nodes
		for (index i = 0; i < numToDel; ++i) {
			// draw a certain (random) number of vertices to be deleted
			node nodeToDel = none;

			do {
				nodeToDel = Aux::Random::integer(G.upperNodeIdBound() - 1);
			} while (! G.hasNode(nodeToDel));

			// mark incident edges first for deletion, delete them outside of neighborhood iterator
			G.forNeighborsOf(nodeToDel, [&](node neigh) {
				edgesToDelete.push_back(std::make_pair(nodeToDel, neigh));
			});
			for (auto elem : edgesToDelete) {
				node& neigh = elem.second;
				G.removeEdge(nodeToDel, neigh);
				GraphEvent event(GraphEvent::EDGE_REMOVAL, nodeToDel, neigh);
//				TRACE("Event: REMOVE edge " , nodeToDel , "-" , neigh);
				eventStream.push_back(event);
			}
			edgesToDelete.clear();

			// eventually delete vertex
			G.removeNode(nodeToDel);
			GraphEvent event(GraphEvent::NODE_REMOVAL, nodeToDel);
//			TRACE("Event: REMOVE node " , nodeToDel);
			eventStream.push_back(event);
		}


		// insert nodes
		for (index i = 0; i < numToIns; ++i) {
			// draw a cluster where the vertex should be inserted, +1 to account for the noise
			count clusterToIns = Aux::Random::integer(
					initGen.numDenseAreas + 1);
			float x; // x-coordinate of new node
			float y; // y-coordinate of new node

			if (clusterToIns < initGen.numDenseAreas) {
				// real cluster, FIXME: DRY!
				// compute random angle between [0, 2pi) and distance between [0, width/2]
				float angle = Aux::Random::real() * 2.0 * PI;
				float dist = Aux::Random::real()
						* initGen.denseAreaXYR[clusterToIns].rad;

				// compute coordinates and adjust them
				x = initGen.denseAreaXYR[clusterToIns].x
						+ cosf(angle) * dist;
				y = initGen.denseAreaXYR[clusterToIns].y
						+ sinf(angle) * dist;
				initGen.moveNodeIntoUnitSquare(x, y);
			} else {
				// noise -> random coordinate
				x = Aux::Random::probability();
				y = Aux::Random::probability();
//				DEBUG("completely random placement");
			}

			// create vertex with these coordinates
			node newNode = G.addNode(x, y);
			Point<float> p(x, y);
			coordinates[newNode] = p;
			GraphEvent event(GraphEvent::NODE_ADDITION, newNode);
//			TRACE("Event: ADD node " , newNode);
			eventStream.push_back(event);
		}


		// determine events by computing new graph structure

		float sqrNeighRad = initGen.neighRad * initGen.neighRad;
		std::map<edge, count> eligibleEdges;

		auto isInRange([&](float squaredDistance) {
			return (squaredDistance <= sqrNeighRad);
		});

		// find for each node the rad-neighborhood
		// FIXME: get rid of quadratic running time!
		G.forNodes([&](node u) {
			std::priority_queue<std::pair<distance, edge> > pq;
			Point<float> p1 = G.getCoordinate(u);
			float& x1 = p1[0];
			float& y1 = p1[1];

			// fill PQ with neighbors in range
			G.forNodes([&](node v) {
				Point<float> p2 = G.getCoordinate(v);
				float& x2 = p2[0];
				float& y2 = p2[1];
				float sqrDist = initGen.squaredDistanceInUnitTorus(x1, y1, x2, y2);

				if (isInRange(sqrDist)) {
					edge e = std::make_pair(std::min(u, v), std::max(u, v));
					pq.push(std::make_pair(-sqrDist, e));
				}
			});

			// mark up to maxNeigh nearest neighbors as eligible
			count end = std::min(initGen.maxNeigh, (count) pq.size());
			for (index i = 0; i < end; ++i) {
				std::pair<distance, edge> currentBest = pq.top();
				pq.pop();
				eligibleEdges[currentBest.second]++;
			}
		});


		// check if edges have to be deleted (was there, but not eligible twice any more)
		G.forEdges([&](node u, node v) {
			edge e = std::make_pair(std::min(u, v), std::max(u, v));
			if (eligibleEdges[e] < 2) {
				edgesToDelete.push_back(std::make_pair(u, v));
			}
			else {
				assert(G.hasEdge(u, v));
				Point<float> p1 = G.getCoordinate(u);
				Point<float> p2 = G.getCoordinate(v);
				edgeweight ew = BASE_WEIGHT / initGen.squaredDistanceInUnitTorus(p1[0], p1[1], p2[0], p2[1]);
				G.setWeight(u, v, ew);
				GraphEvent event(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, ew);
//				TRACE("Event: UPD edge weight " , u , "-" , v, ", weight: ", ew);
				eventStream.push_back(event);
			}
			eligibleEdges.erase(e);
		});

		// now really delete the edges after leaving the G.forEdges iterator
		for (auto elem : edgesToDelete) {
			node& u = elem.first;
			node& v = elem.second;
			G.removeEdge(u, v);
			GraphEvent event(GraphEvent::EDGE_REMOVAL, u, v);
//			TRACE("Event: REMOVE edge " , u , "-" , v);
			eventStream.push_back(event);
		}
		edgesToDelete.clear();


		// check if edges have to be inserted
		for (auto edgePair : eligibleEdges) {
			if (edgePair.second >= 2) {
				node u = edgePair.first.first;
				node v = edgePair.first.second;
				assert(! G.hasEdge(u, v) && G.hasNode(u) && G.hasNode(v));

				Point<float> p1 = G.getCoordinate(u);
				Point<float> p2 = G.getCoordinate(v);
				edgeweight ew = BASE_WEIGHT / initGen.squaredDistanceInUnitTorus(p1[0], p1[1], p2[0], p2[1]);
				G.addEdge(u, v, ew);
				assert(G.hasEdge(u, v));
				GraphEvent event(GraphEvent::EDGE_ADDITION, u, v, ew);
//				TRACE("Event: ADD edge " , u , "-" , v, ", weight: ", ew);
				eventStream.push_back(event);
			}
		}


		eventStream.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}

//	DGSWriter dgsWriter;
//	dgsWriter.write(eventStream, "output/eventStream.dgs");
//	for (auto event : eventStream) {
//		TRACE("event: ", event.type, ", node: ", event.u);
//	}

	return eventStream;
}

std::map<node, Point<float> > DynamicPubWebGenerator::getNewCoordinates() const {
	return coordinates;
}

} /* namespace NetworKit */
