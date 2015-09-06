/*
 * KruskalMSF.h
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#ifndef KRUSKALMSF_H_
#define KRUSKALMSF_H_

#include "../Globals.h"
#include "Graph.h"
#include "SpanningForest.h"

namespace NetworKit {

class KruskalMSF: public SpanningForest {
public:
	KruskalMSF(const Graph& G);
	virtual ~KruskalMSF() = default;

	virtual void run() override;
};

} /* namespace NetworKit */
#endif /* KRUSKALMSF_H_ */
