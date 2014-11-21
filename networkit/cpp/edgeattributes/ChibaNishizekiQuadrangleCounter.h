/*
 * ChibaNishizekiQuadrangleCounter.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef CHIBANISHIZEKIQUADRANGLECOUNTER_H
#define CHIBANISHIZEKIQUADRANGLECOUNTER_H

#include "EdgeAttribute.h"

namespace NetworKit {

class ChibaNishizekiQuadrangleCounter : public EdgeAttribute<int> {

protected:
	const Graph& G;

public:
	ChibaNishizekiQuadrangleCounter(const Graph& G);

	virtual std::vector<int> getAttribute() override;
};

} // namespace NetworKit

#endif // CHIBANISHIZEKIQUADRANGLECOUNTER_H
