/*
 * ChibaNishizekiTriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef CHIBANISHIZEKI_H_
#define CHIBANISHIZEKI_H_

#include "../graph/Graph.h"
#include "EdgeAttribute.h"

namespace NetworKit {

/**
 * An implementation of the triangle counting algorithm by Chiba/Nishizeki.
 */
class ChibaNishizekiTriangleCounter : public EdgeAttribute<count> {

protected:
	const Graph& G;

public:

	ChibaNishizekiTriangleCounter(const Graph& G);

	virtual std::vector<count> getAttribute() override;
};

} /* namespace NetworKit */

#endif /* CHIBANISHIZEKI_H_ */
