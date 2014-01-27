/*
 * MaxentStress.cpp
 *
 *  Created on: 22.01.2014
 *      Author: Henning
 */

#include "MaxentStress.h"

namespace NetworKit {

MaxentStress::MaxentStress(Point<float> bottomLeft, Point<float> topRight, bool useGivenLayout):
				Layouter(bottomLeft, topRight, useGivenLayout)
{

}

MaxentStress::~MaxentStress() {

}

void MaxentStress::draw(Graph& G) {
	count n = G.numberOfNodes();
	initialize(G);

	//////////////////////////////////////////////////////////
	// Force calculations
	//////////////////////////////////////////////////////////



	//////////////////////////////////////////////////////////
	// Move vertices according to forces
	//////////////////////////////////////////////////////////
	auto move([&](Point<float>& p, Point<float>& force, float step) {
		// x_i := x_i + step * (f / ||f||)
		float len = force.length();
		if (len > 0) {
			p += force.scale(step / len);
		}

		// position inside frame
		p[0] = fmax(p[0], 0.0);
		p[1] = fmax(p[1], 0.0);
		p[0] = fmin(p[0], 1.0);
		p[1] = fmin(p[1], 1.0);
	});


	//////////////////////////////////////////////////////////
	// Cooling schedule
	//////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////
	// Check convergence
	//////////////////////////////////////////////////////////
	auto isConverged([&](std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		float change = 0.0;

		for (index i = 0; i < oldLayout.size(); ++i) {
			change += oldLayout[i].distance(newLayout[i]); // could be accelerated by squared distance
		}
		DEBUG("change: ", change);

		return (change < 0.1); // FIXME: externalize
	});



	//////////////////////////////////////////////////////////
	// Preparations for main loop
	//////////////////////////////////////////////////////////
	bool converged = false;
	std::vector<float> origin = {0.0, 0.0};
	std::vector<Point<float> > forces(n, origin);
	count iter = 0;

	//////////////////////////////////////////////////////////
	// Main loop
	//////////////////////////////////////////////////////////

	// Paper: q = 0 for many graphs, q = 0.8 for graphs with many degree-1 nodes
	// alpha: initially 1, then in each iteration alpha := 0.3 * alpha
	float q = 0;
	float alpha = 1.0; // FIXME: make alpha dependent on iteration

	while (! converged) {
		std::vector<Point<float> > previousLayout = layout;

		// init for current iteration
		G.parallelForNodes([&](node u) {
			forces[u] = origin;
		});

		// apply forces to each node
		G.forNodes([&](node u) {
			// FIXME: take care of singletons...
			assert(G.weightedDegree(u) > 0.0);

			Point<float> uPoint = layout[u];
			Point<float> attractiveForce(0.0, 0.0);
			Point<float> repulsiveForce(0.0, 0.0);

			G.forNodes([&](node v) {
				if (u < v) { // only unordered pairs
					Point<float> vPoint = layout[v];
					float diffX = uPoint[0] - vPoint[0];
					float diffY = uPoint[1] - vPoint[1];
					Point<float> diffVec(diffX, diffY);
					float len = diffVec.length();

					if (len > 0.0) {
						if (G.hasEdge(u, v)) {
							// sum over all node pairs in S
							// $\frac{1}{\rho_i} \sum_{i,j \in S} w_{ij} * (x_j + d_{ij} \frac{x_i - x_j}{\Vert x_i - x_j \Vert})$
								float dist = 0.1; // FIXME: dist = ... (algebraic distance between u and v)
								diffVec.scale(dist / len);
								diffVec += vPoint;
								attractiveForce += diffVec.scale(G.weight(u, v));
						}
						else {
							// traverse remaining vertices not in neighborhood for repulsive forces
							// sum over all node pairs not in S
							// $\frac{\alpha}{\rho_i} \sum_{i,j \notin S} w_{ij} \frac{x_i - x_j}{\Vert x_i - x_j \Vert^{q+2}}$
							// TODO: approximation (e.g. Barnes-Hut)

							float denom = 1.0 / pow(len, q+2);
							diffVec.scale(denom);
							repulsiveForce += diffVec.scale(G.weight(u, v));
						}
					}
				}
			});

			// apply forces to node u
			float rhoInv = 1.0 / (float) G.weightedDegree(u);
			attractiveForce.scale(rhoInv);
			repulsiveForce.scale(alpha * rhoInv);

			forces[u] += attractiveForce;
			forces[u] += repulsiveForce;
		});


		// move nodes
		G.forNodes([&](node u) { // FIXME: parallel!
			move(layout[u], forces[u], 1.0); // FIXME: step length

			TRACE("moved ", u, " by: ", forces[u][0], " and ", forces[u][1]);
			TRACE("old pos: ", previousLayout[u].toString(), ", new pos: ", layout[u].toString());
		});

		++iter;
		converged = isConverged(previousLayout, layout) || iter >= 500; // FIXME: externalize
		DEBUG("iteration finished: ", iter, "converged: ", converged);
	}

	// copy layout into graph
	G.parallelForNodes([&](node u) {
		for (index d = 0; d < layout[u].getDimensions(); ++d) { // TODO: accelerate loop, needs Point<T> as coordinate type in graph.h
			G.setCoordinate(u, layout[u]);
		}
	});
}

} /* namespace NetworKit */
