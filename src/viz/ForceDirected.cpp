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
	// TODO Auto-generated constructor stub

	// FIXME: set forceCorrection in base class
	forceCorrection = 0.05;
}

ForceDirected::~ForceDirected() {
	// TODO Auto-generated destructor stub
}

void ForceDirected::draw(Graph& g) {
	const float INITIAL_STEP_LENGTH = 1.0;

	int width = (topRight.getValue(0) - bottomLeft.getValue(0));
	int height = (topRight.getValue(1) - bottomLeft.getValue(1));
	int area = width * height;
	count n = g.numberOfNodes();

	// initialize randomly
	randomInitCoordinates(g);


	auto attractiveForce([&](Point<float>& p1, Point<float>& p2) {
		Point<float> force = p1;

		float dist = p1.distance(p2);
		float strength = dist / forceCorrection;
		force.scale(strength); // FIXME:

		return force;
	});

	auto repellingForce([&](Point<float>& p1, Point<float>& p2) {
		Point<float> force = p1;

		float sqrDist = p1.squaredDistance(p2);
		float strength = forceCorrection * forceCorrection / sqrDist; // TODO: store square
		force.scale(strength);

		return force;
	});

	auto move([&](Point<float>& p, Point<float>& force, float step) {
		Point<float> newPoint = p;

		// x_i := x_i + step * (f / ||f||)
		p += force.scale(step / force.length());

		return newPoint;
	});

	auto isConverged([&](std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		float eps = 0.1;
		float change = 0.0;

		for (index i = 0; i < oldLayout.size(); ++i) {
			change += oldLayout[i].distance(newLayout[i]);
		}

		TRACE("change: " << change);

		return (change < eps);
	});

	auto updateStepLength([&](float step, std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		float newStep = 0.95;
		// FIXME: clever update
		return newStep;
	});

	float optimalPairwiseDist = sqrt((float) area / (float) n);
	bool converged = false;
	float step = INITIAL_STEP_LENGTH;
	std::vector<float> origin = {0.0, 0.0};
	count numIter = 0;

	converged = true; // FIXME: just for test
	while (! converged) {
		std::vector<Point<float> > previousLayout = layout;
		g.forNodes([&](node u) {
			Point<float> force = origin;

			g.forNeighborsOf(u, [&](node v) {
				force += attractiveForce(layout[u], layout[v]);
			});

			g.forNodes([&](node v) {
				if (u != v) {
					force += repellingForce(layout[u], layout[v]);
				}
			});

			move(layout[u], force, step);
		});

		++numIter;
		step = updateStepLength(step, previousLayout, layout);
		converged = isConverged(previousLayout, layout) || numIter >= MAX_ITER;
	}

	// copy layout into graph
	g.forNodes([&](node u) {
		for (index d = 0; d < layout[u].getDimensions(); ++d) { // TODO: accelerate loop
			g.setCoordinate(u, d, layout[u].getValue(d));
		}
	});
}



} /* namespace NetworKit */
