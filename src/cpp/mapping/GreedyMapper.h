/*
 * GreedyMapper.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef GREEDYMAPPER_H_
#define GREEDYMAPPER_H_

#include "StaticMapper.h"

namespace NetworKit {

/**
 * TODO: class documentation
 */
class GreedyMapper: public NetworKit::StaticMapper {
public:
	GreedyMapper();
	virtual ~GreedyMapper();

	//virtual Mapping run(Graph& guest, Graph& host);
};

} /* namespace NetworKit */
#endif /* GREEDYMAPPER_H_ */
