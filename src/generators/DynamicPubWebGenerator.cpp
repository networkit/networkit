/*
 * DynamicPubWebGenerator.cpp
 *
 *  Created on: 15.01.2014
 *      Author: Henning
 */

#include "DynamicPubWebGenerator.h"

namespace NetworKit {

DynamicPubWebGenerator::DynamicPubWebGenerator(count numNodes,
		count numberOfDenseAreas, float neighborhoodRadius,
		count maxNumberOfNeighbors) :
		initGen(numNodes, numberOfDenseAreas, neighborhoodRadius,
				maxNumberOfNeighbors) {
	G = initGen.generate();
}

DynamicPubWebGenerator::~DynamicPubWebGenerator() {

}

std::vector<GraphEvent> DynamicPubWebGenerator::generate(count nSteps) {

	count numToDel = (count) (G.numberOfNodes() * 0.05); // TODO: externalize, possibly randomize
	count numToIns = (count) (G.numberOfNodes() * 0.05); // TODO: externalize, possibly randomize
	std::vector<GraphEvent> eventStream;
	coordinates.clear();

	for (index step = 0; step < nSteps; ++step) {

		// delete nodes
		for (index i = 0; i < numToDel; ++i) {
			// draw a certain (random) number of vertices to be deleted
			node nodeToDel = Aux::Random::integer(G.upperNodeIdBound() - 1);

			if (G.hasNode(nodeToDel)) {
				// delete incident edges first
				G.forNeighborsOf(nodeToDel, [&](node neigh) {
					G.removeEdge(nodeToDel, neigh);
					GraphEvent event(GraphEvent::EDGE_REMOVAL, nodeToDel, neigh);
					TRACE("Event: REMOVE edge " << nodeToDel << "-" << neigh);
					eventStream.push_back(event);
				});

				// eventually delete vertex
				G.removeNode(nodeToDel);
				GraphEvent event(GraphEvent::NODE_REMOVAL, nodeToDel);
				TRACE("Event: REMOVE node " << nodeToDel);
				eventStream.push_back(event);
			}
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
				float angle = Aux::Random::real() * 2.0 * M_PI;
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
			}

			// create vertex with these coordinates
			node newNode = G.addNode(x, y);
			Point<float> p(x, y);
			coordinates[newNode] = p;
			GraphEvent event(GraphEvent::NODE_ADDITION, newNode);
			TRACE("Event: ADD node " << newNode);
			eventStream.push_back(event);
		}

		// determine events by computing new graph structure

		float sqrNeighRad = initGen.neighRad * initGen.neighRad;
		std::map<edge, count> eligibleEdges;

		auto isInRange([&](float squaredDistance) {
			return (squaredDistance <= sqrNeighRad);
		});

		G.forNodes([&](node u) {
			std::priority_queue<std::pair<distance, edge> > pq;
			float x1 = G.getCoordinate(u, 0);
			float y1 = G.getCoordinate(u, 1);

			// fill PQ with neighbors in range
			G.forNodes([&](node v) {
				float x2 = G.getCoordinate(v, 0);
				float y2 = G.getCoordinate(v, 1);
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
				G.removeEdge(u, v);
				GraphEvent event(GraphEvent::EDGE_REMOVAL, u, v);
				TRACE("Event: REMOVE edge " << u << "-" << v);
				eventStream.push_back(event);
			}
			eligibleEdges.erase(e);
		});

		// check if edges have to be inserted
		for (auto edgePair : eligibleEdges) {
			if (edgePair.second >= 2) {
				node u = edgePair.first.first;
				node v = edgePair.first.second;
				// TODO: make common routine
				edgeweight ew = BASE_WEIGHT / initGen.squaredDistanceInUnitTorus(
						G.getCoordinate(u, 0), G.getCoordinate(u, 1), G.getCoordinate(v, 0), G.getCoordinate(v, 1));
				G.addEdge(u, v);
				GraphEvent event(GraphEvent::EDGE_ADDITION, u, v, ew);
				TRACE("Event: ADD edge " << u << "-" << v);
				eventStream.push_back(event);
			}
		}

		eventStream.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}

	// TODO: remove this test
	DGSWriter dgsWriter;
	dgsWriter.write(eventStream, "output/eventStream.dgs");

	return eventStream;
}

std::map<node, Point<float> > DynamicPubWebGenerator::getNewCoordinates() const {
	return coordinates;
}

} /* namespace NetworKit */
