/*
 * Overlapper.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef OVERLAPPER_H_
#define OVERLAPPER_H_

#include <set>
#include <vector>

#include "../graph/Graph.h"
#include "../clustering/base/Clustering.h"

namespace EnsembleClustering {


class Overlapper {

public:

	Overlapper();

	virtual ~Overlapper();

	virtual Clustering run(Graph& G, std::vector<Clustering>& clusterings) = 0;
};

} /* namespace EnsembleClustering */
#endif /* OVERLAPPER_H_ */
