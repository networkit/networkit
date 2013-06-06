/*
 * ForceDirected.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#include "ForceDirected.h"

namespace NetworKit {

ForceDirected::ForceDirected(Point<float> bottom_left, Point<float> top_right):
		SpringEmbedder(bottom_left, top_right) {

}

ForceDirected::~ForceDirected() {

}

void ForceDirected::draw(Graph& g) {
	const float INITIAL_STEP_LENGTH = 1.0;

	int width = (topRight.getValue(0) - bottomLeft.getValue(0));
	int height = (topRight.getValue(1) - bottomLeft.getValue(1));
	int area = width * height;
	count n = g.numberOfNodes();
	double optPairDist = 0.95 * sqrt((double) area / (double) n);
	TRACE("k: " << optPairDist);

	// initialize randomly
	randomInitCoordinates(g);


	auto attractiveForce([&](Point<float>& p1, Point<float>& p2) {
		Point<float> force = p1 - p2;

		float dist = force.length();
		float strength = dist / optPairDist;
		force.scale(-strength);

		return force;
	});



	auto repellingForce([&](Point<float>& p1, Point<float>& p2) {
		Point<float> force = p1 - p2;

		float sqrDist = force.squaredLength();
		float strength = optPairDist * optPairDist / sqrDist; // TODO: store square
		force.scale(strength);

		return force;
	});

	auto move([&](Point<float>& p, Point<float>& force, float step) {
		// x_i := x_i + step * (f / ||f||)
		p += force.scale(step / force.length());

		// position inside frame
		p.setValue(0, fmax(p.getValue(0), 0.0));
		p.setValue(1, fmax(p.getValue(1), 0.0));
		p.setValue(0, fmin(p.getValue(0), 1.0));
		p.setValue(1, fmin(p.getValue(1), 1.0));
	});

	auto isConverged([&](std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		float eps = 1e-1;
		float change = 0.0;

		for (index i = 0; i < oldLayout.size(); ++i) {
			change += oldLayout[i].distance(newLayout[i]);
		}

		TRACE("change: " << change);

		return (change < eps);
	});

	auto updateStepLength([&](float step, std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		float newStep = 1.0 / step;
		newStep += 0.5;

		return 1.0 / newStep;
	});


	bool converged = false;
	float step = INITIAL_STEP_LENGTH;
	std::vector<float> origin = {0.0, 0.0};
	count numIter = 0;
	std::vector<Point<float> > forces(n, origin);

	while (! converged) {
		std::vector<Point<float> > previousLayout = layout;

		// repulsive forces
		g.forNodes([&](node u) {
			forces[u] = origin;
			g.forNodes([&](node v) {
				if (u != v) {
					forces[u] += repellingForce(previousLayout[u], previousLayout[v]);
				}
			});
		});


		// attractive forces
		g.forEdges([&](node u, node v) {
			// TOOD: store result instead of computing it twice; requires point operations
			forces[u] += attractiveForce(previousLayout[u], previousLayout[v]);
			forces[v] += attractiveForce(previousLayout[v], previousLayout[u]);
		});


		// move nodes
		g.forNodes([&](node u) {
			move(layout[u], forces[u], step);

//			DEBUG("moved " << u);
//			DEBUG("by: " << forces[u].getValue(0) << " and " << forces[u].getValue(1));
//			DEBUG("new x: " << layout[u].getValue(0));
//			DEBUG("new y: " << layout[u].getValue(1));
		});

		++numIter;
		step = updateStepLength(step, previousLayout, layout);
		converged = isConverged(previousLayout, layout) || numIter >= MAX_ITER;

		TRACE("Iteration finished: " << numIter);
	}

	// copy layout into graph
	g.forNodes([&](node u) {
		for (index d = 0; d < layout[u].getDimensions(); ++d) { // TODO: accelerate loop
			g.setCoordinate(u, d, layout[u].getValue(d));
		}
	});
}



} /* namespace NetworKit */
