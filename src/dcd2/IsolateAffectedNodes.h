/*
 * IsolateAffectedNodes.h
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#ifndef ISOLATEAFFECTEDNODES_H_
#define ISOLATEAFFECTEDNODES_H_

#include "PrepStrategy2.h"

namespace NetworKit {

class IsolateAffectedNodes: public PrepStrategy2 {

	virtual void process(Partition& communities, std::vector<GraphEvent> eventStream);

};

} /* namespace NetworKit */

#endif /* ISOLATEAFFECTEDNODES_H_ */
