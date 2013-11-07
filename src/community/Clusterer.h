/*
 * Clusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERER_H_
#define CLUSTERER_H_

#include "../clustering/Clustering.h"

namespace NetworKit {


class Clusterer {
public:

	Clusterer();

	virtual ~Clusterer();

	virtual Clustering run(Graph& G) = 0;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif /* CLUSTERER_H_ */
