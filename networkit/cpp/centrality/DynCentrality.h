/*
 * DynSSSP.h
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

#ifndef DYNCENTRALITY_H_
#define DYNCENTRALITY_H_

#include "../dynamics/GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Interface for dynamic centrality algorithms.
 */
class DynCentrality {

public:

    virtual void update(const std::vector<GraphEvent>& batch) = 0;

};

} /* namespace NetworKit */

#endif /* DYNCENTRALITY_H_ */
