/*
 * DynamicHyperbolicGenerator.h
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#ifndef DYNAMICHYPERBOLICGENERATOR_H_
#define DYNAMICHYPERBOLICGENERATOR_H_

#include <map>

#include "DynamicGraphGenerator.h"
#include "quadtree/Quadtree.h"


namespace NetworKit {

class DynamicHyperbolicGenerator: public NetworKit::DynamicGraphGenerator  {
	friend class GeneratorsGTest;
public:
	/**
	 * Initialize a dynamic hyperbolic generator and generate initial node positions
	 * Node movement happens in discrete time steps
	 *
	 * @param n number of nodes
	 * @param avgDegree expected average degree of target graph
	 * @param exp exponent of power-law degree distribution
	 * @param T temperature parameter in edge probabilities
	 * @param moveEachStep fraction of nodes which are moved at each time step, should be non-negative
	 * @param moveDistance base value for the node movements
	 */

	DynamicHyperbolicGenerator(count n = 1000, double avgDegree=6, double exp=3, double T=0, double moveEachStep = 0, double moveDistance = 0);

	/**
	 * Initialize a dynamic hyperbolic generator with given initial node positions in polar coordinates
	 * Node movement happens in discrete time steps
	 *
	 * @param angles angular coordinates of initial positions
	 * @param radii radial coordinates of initial positions
	 * @param R radius of hyperbolic disk
	 * @param alpha dispersion parameter of point distribution
	 * @param T temperature parameter in edge probabilities
	 * @param moveEachStep fraction of nodes which are moved at each time step, should be non-negative
	 * @param moveDistance base value for the node movements
	 */
	DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii,  double R, double alpha, double T=0, double moveEachStep = 0, double moveDistance = 0);

	/**
	 * Default constructor
	 */
	DynamicHyperbolicGenerator();

	/**
	 * Run the dynamic, changes state of generator
	 *
	 * @param nSteps number of time steps to iterate over
	 * @return changed edges
	 */
	std::vector<GraphEvent> generate(count nSteps) override;

	/**
	 * Get the graph corresponding to the current state of the generator. Does not change the generator
	 *
	 * @return graph at the current state
	 */
	Graph getGraph() const;

	/**
	 * Get coordinates in native representation
	 * @return vector of 2D-Points in Cartesian coordinates
	 */
	std::vector<Point<float> > getCoordinates() const;

private:
	/**
	 * Generate initial node positions and fill the quadtree with them
	 */
	void initializePoints();
	void initializeQuadTree();

	/**
	 * Generate initial movement vectors for all points
	 */
	void initializeMovement();

	/**
	 * Generate initial movement vectors for all points
	 */
	void recomputeBands();

	vector<index> getNeighborsInBands(index i, bool bothDirections=true);

	/**
	 * Execute node movement part of time step
	 *
	 * @param result vector to store GraphEvents in
	 */
	void getEventsFromNodeMovement(vector<GraphEvent> &result);

	/**
     * Move a single node
	 *
	 * @param node Index of the node that should be moved.
	 */
	void moveNode(index node);

	//general geometry parameters
	count nodeCount;
	double alpha;
	double R;
	double T;

	//movement parameters
	double moveEachStep;
	double moveDistance;

	//coordinates
	vector<double> angles;
	vector<double> radii;

	//movement vectors
	vector<double> angularMovement;
	vector<double> radialMovement;

	//data structures
	Quadtree<index, false> quad;
	vector<double> bandRadii;
	vector<vector<Point2D<double>>> bands;
	vector<vector<double> > bandAngles;

	bool initialized;
};

} /* namespace NetworKit */
#endif /* DYNAMICHYPERBOLICGENERATOR_H_ */

