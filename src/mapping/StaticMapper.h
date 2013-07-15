/*
 * StaticMapper.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef STATICMAPPER_H_
#define STATICMAPPER_H_

#include <map>
#include "../graph/Graph.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

class StaticMapper {
public:
	StaticMapper();
	virtual ~StaticMapper();

	virtual std::map<index, index> run(Graph& guest, Graph& host) = 0;
};

} /* namespace NetworKit */
#endif /* STATICMAPPER_H_ */
