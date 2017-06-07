#ifndef DYNALGORITHM_H
#define DYNALGORITHM_H

#include <string>
#include <stdexcept>
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

class DynAlgorithm {

public:
	/**
	 * Virtual default destructor
	 */
	virtual ~DynAlgorithm() = default;

	/**
	 * The generic update method for updating data structure after an update.
	 */
	virtual void update(GraphEvent e) = 0;

	/**
	 * The generic update method for updating data structure after a batch of updates.
	 */
	virtual void updateBatch(const std::vector<GraphEvent>& batch) = 0;

};

} /* NetworKit */

#endif /* DYNALGORITHM_H */
