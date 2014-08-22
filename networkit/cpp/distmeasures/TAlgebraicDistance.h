/*
 * TAlgebraicDistance.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TALGEBRAICDISTANCE_H_
#define TALGEBRAICDISTANCE_H_

#include "TNodeDistance.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distmeasures
 * Algebraic distance assigns a distance value to pairs of nodes
 * according to their structural closeness in the graph. 
 */
class TAlgebraicDistance: public NetworKit::TNodeDistance {

public:

	TAlgebraicDistance(const Graph& G);

	/** Default destructor */
	virtual ~TAlgebraicDistance();

	void initialize(const Parameters& param);

	double distance(node u, node v);

protected:

	count numSystems; //!< number of vectors/systems used for algebraic iteration
	count numIters; //!< number of iterations in each system
	double omega; //!<
	index norm;
	const index MAX_NORM = 0;

	std::vector<std::vector<double> > loads; //!< loads[i]: vector of loads of length n for one system


	void randomInit();


};

} /* namespace NetworKit */
#endif /* TALGEBRAICDISTANCE_H_ */
