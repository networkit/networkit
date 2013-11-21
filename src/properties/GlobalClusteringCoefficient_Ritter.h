/*
 * GlobalClusteringCoefficient_Ritter.h
 */

#ifndef GLOBALCLUSTERINGCOEFFICIENT_RITTER_H_
#define GLOBALCLUSTERINGCOEFFICIENT_RITTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class GlobalClusteringCoefficient_Ritter {

public:

	GlobalClusteringCoefficient_Ritter();

	virtual ~GlobalClusteringCoefficient_Ritter();

	virtual double calculate(Graph& G);
};

}
#endif /* GLOBALCLUSTERINGCOEFFICIENT_RITTER_H_ */


