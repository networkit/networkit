/*
 * PrepStrategy2.h
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#ifndef PREPSTRATEGY2_H_
#define PREPSTRATEGY2_H_

#include "../structures/Partition.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

class PrepStrategy2 {

public:

	virtual void process(Partition& communities, std::vector<GraphEvent> eventStream) = 0;

};

} /* namespace NetworKit */

#endif /* PREPSTRATEGY2_H_ */
