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
#include "../graph/GraphDistance.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

typedef std::vector<index> Permutation;
typedef std::map<index, index> Mapping;

class StaticMapper {
protected:

public:
	StaticMapper();
	virtual ~StaticMapper();

	virtual Mapping run(Graph& guest, Graph& host) = 0;

	virtual Mapping trivial(Graph& guest, Graph& host);

	virtual edgeweight cost(const Graph& guest, const Graph& host, Mapping& mapping);
};

} /* namespace NetworKit */
#endif /* STATICMAPPER_H_ */
