/*
 * ClusteringCoefficient.h
 *
 *  Created on: 08.04.2013
 *      Author: cls
 */

#ifndef GLOBALCLUSTERINGCOEFFICIENT_H_
#define GLOBALCLUSTERINGCOEFFICIENT_H_

#include "../graph/Graph.h"

namespace NetworKit {

// TODO: is class necessary?
class GlobalClusteringCoefficient {

public:

	GlobalClusteringCoefficient();

	virtual ~GlobalClusteringCoefficient();

	virtual float run(Graph& G);
};

}
#endif




























