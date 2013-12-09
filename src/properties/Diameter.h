/*
 * Diameter.h
 *
 *  Created on: 21.11.13
 *      Author: lbarth, dweiss
 */

#ifndef DIAMETER_H_
#define DIAMETER_H_

#include "../graph/Graph.h"

namespace GrauBart {

using namespace NetworKit;

// TODO: is class necessary?
class Diameter {

private:
  count getUpperBound(const Graph &G, const node root) const;

public:

	Diameter();

	virtual ~Diameter();

	virtual std::pair<count,count> estimateDiameterRange(const Graph& G) const;  
};

} /* namespace Graubart */
#endif /* DIAMETER_H_ */
