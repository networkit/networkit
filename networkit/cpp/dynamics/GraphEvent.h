/*
 * GraphEvent.h
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENT_H_
#define GRAPHEVENT_H_

#include "../graph/Graph.h"

namespace NetworKit {


/**
 * @ingroup dynamics
 */
class GraphEvent {


public:

	enum Type {
		NODE_ADDITION,
		NODE_REMOVAL,
		NODE_RESTORATION,
		EDGE_ADDITION,
		EDGE_REMOVAL,
		EDGE_WEIGHT_UPDATE, 
		EDGE_WEIGHT_INCREMENT,
		TIME_STEP
		
	};

	Type type;	//!< type of graph event
	node u; 				//!< first node parameter
	node v;					//!< second node parameter
	edgeweight w;			//!< edge weight parameter


	GraphEvent() = default;

	GraphEvent(Type type, node u = none, node v = none, edgeweight w = 1.0);

	static bool compare(GraphEvent a, GraphEvent b);
	static bool equal(GraphEvent a, GraphEvent b);

	/**
	 * Return string representation.
	 */
	std::string toString();

};

} /* namespace NetworKit */
#endif /* GRAPHEVENT_H_ */
