/*
 * RcmMapper.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef RCMMAPPER_H_
#define RCMMAPPER_H_

#include "StaticMapper.h"

namespace NetworKit {

typedef std::vector<index> Permutation;
typedef std::map<index, index> Mapping;

class RcmMapper: public NetworKit::StaticMapper {
public:
	RcmMapper();
	virtual ~RcmMapper();

	virtual Mapping run(Graph& guest, Graph& host);

	Permutation permute(const Graph& graph) const;
	Permutation invert(const Permutation& piIn) const;
};

} /* namespace NetworKit */
#endif /* RCMMAPPER_H_ */
