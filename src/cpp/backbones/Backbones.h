/*
 * Backbones.h
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef BACKBONES_H_
#define BACKBONES_H_

#include "../graph/Graph.h"
#include "AttributeGenerator.h"
#include "ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

/** 
 * Combines attribute generators and edge attribute filters into different
 * backbone algorithms.
 */

class SimmelianBackboneParametric : public BackboneCalculator {
	Graph SimmelianBackboneParametric::calculate(const Graph& g, const edgeAttribute& attribute) {
		ChibaNishizekiTriangleCounter counter;
		edgeAttribute triangles = counter.getAttribute(g);

		//TODO: ...

		return g;
	}
};


} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
