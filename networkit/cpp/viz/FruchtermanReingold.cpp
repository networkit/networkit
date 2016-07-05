/*
 * ForceDirected.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#include "FruchtermanReingold.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

const float FruchtermanReingold::INITIAL_STEP_LENGTH = 1.0;
const float FruchtermanReingold::OPT_PAIR_SQR_DIST_SCALE = 0.3;


FruchtermanReingold::FruchtermanReingold(Point<float> bottom_left, Point<float> top_right, bool useGivenCoordinates, count maxIterations, float precision):
		Layouter(bottom_left, top_right, useGivenCoordinates), maxIter(maxIterations), prec(precision), step(INITIAL_STEP_LENGTH)
{

}

void FruchtermanReingold::draw(Graph& g) {

	int width = (topRight[0] - bottomLeft[0]);
	int height = (topRight[1] - bottomLeft[1]);
	int area = width * height;
	count n = g.numberOfNodes();

	float optPairSqrDist = OPT_PAIR_SQR_DIST_SCALE * (float) area / (float) n;
	float optPairDist = sqrt(optPairSqrDist);
	DEBUG("optPairDist: ", optPairDist);

	initialize(g);

	//////////////////////////////////////////////////////////
	// Force calculations
	//////////////////////////////////////////////////////////
	auto attractiveForce([&](Point<float>& p1, Point<float>& p2) {
		Point<float> force = p1 - p2;

		float dist = force.length();
		float strength = dist / optPairDist;
		force.scale(strength);

		return force;
	});

	auto repulsiveForce([&](Point<float>& p1, Point<float>& p2) {
		Point<float> force = p1 - p2;

		float sqrDist = force.squaredLength();
		if (sqrDist > 0) {
			float strength = optPairSqrDist / sqrDist;
			force.scale(strength);
		}

		return force;
	});



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
	auto updateStepLength([&](std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		step += 0.1; // TODO: externalize
		return 1.0 / step;
	});


	//////////////////////////////////////////////////////////
	// Check convergence
	//////////////////////////////////////////////////////////
	auto isConverged([&](std::vector<Point<float> >& oldLayout,
			std::vector<Point<float> >& newLayout) {
		float change = g.parallelSumForNodes([&](node i) {
			return oldLayout[i].distance(newLayout[i]); // could be accelerated by squared distance
		});
		DEBUG("change: ", change);

		return (change < prec);
	});



	//////////////////////////////////////////////////////////
	// Preparations for main loop
	//////////////////////////////////////////////////////////
	bool converged = false;
	std::vector<float> origin = {0.0, 0.0};
	std::vector<Point<float> > forces(n, origin);
	float actualStep = INITIAL_STEP_LENGTH;
	count iter = 0;

	//////////////////////////////////////////////////////////
	// Main loop
	//////////////////////////////////////////////////////////
	while (! converged) {
		std::vector<Point<float> > previousLayout = layout;

		// init for current iteration
		g.parallelForNodes([&](node u) {
			forces[u] = origin;
		});

		// repulsive forces
		g.parallelForNodePairs([&](node u, node v) {
			Point<float> force = repulsiveForce(previousLayout[u], previousLayout[v]);
			forces[u] += force;
			forces[v] -= force;
		});

		// attractive forces
		g.parallelForEdges([&](node u, node v) {
			Point<float> attr = attractiveForce(previousLayout[u], previousLayout[v]);
			forces[u] -= attr;
			forces[v] += attr;
		});


		// move nodes
		g.parallelForNodes([&](node u) {
			move(layout[u], forces[u], actualStep);

//			TRACE("moved ", u, " by: ", forces[u][0], " and ", forces[u][1]);
//			TRACE("old pos: ", previousLayout[u].toString(), ", new pos: ", layout[u].toString());
		});

		++iter;
		actualStep = updateStepLength(previousLayout, layout);
		converged = isConverged(previousLayout, layout) || iter >= maxIter;

		DEBUG("new step length: ", actualStep, ", iteration finished: ", iter);
	}

	// copy layout into graph
	g.parallelForNodes([&](node u) {
		g.setCoordinate(u, layout[u]);
		DEBUG("coordinate of ", u, ": ", layout[u].toString());
	});
}



} /* namespace NetworKit */
