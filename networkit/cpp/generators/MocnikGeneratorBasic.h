/*
 * MocnikGeneratorBasic.h
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#ifndef MOCNIKGENERATORBASIC_H_
#define MOCNIKGENERATORBASIC_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

/**
 * @ingroup generators
 */
class MocnikGeneratorBasic: public StaticGraphGenerator {
private:
	// GENERAL DATA
	
	/**
	 * Position of each node in space.  The index of the vector is also the number of
	 * the node.
	 */
	std::vector<std::vector<double>> nodePositions;

protected:
	count dim;
	count n;
	double k;

public:
	/**
	 * Creates random spatial graphs according to the Mocnik model.
	 *
	 * Please cite the following publications, in which you will find a
	 * description of the model:
	 *
	 * Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
	 * the Context of Local and Global Optimization", Scientific Reports 8(11274)
	 * 2018. doi: 10.1038/s41598-018-29131-0
	 *
	 * Franz-Benjamin Mocnik, Andrew Frank: "Modelling Spatial Structures",
	 * Proceedings of the 12th Conference on Spatial Information Theory (COSIT),
	 * 2015, pages 44-64. doi: 10.1007/978-3-319-23374-1_3
	 *
	 * Non-improved algorithm.
	 *
	 * @param dim  Dimension of the space.
	 * @param n  Number of nodes in the graph.
	 * @param k  Density parameter, determining the ratio of edges to nodes.
	 */
	MocnikGeneratorBasic(count dim, count n, double k);

	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* MOCNIKGENERATORBASIC_H_ */
