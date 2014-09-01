/*
 * GeneratorsGTest.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef GENERATORSGTEST_H_
#define GENERATORSGTEST_H_

#include <gtest/gtest.h>
#include <cmath>

#include "../DynamicGraphSource.h"
#include "../DynamicBarabasiAlbertGenerator.h"
#include "../PubWebGenerator.h"
#include "../DynamicPubWebGenerator.h"
#include "../HyperbolicGenerator.h"
#include "../DynamicHyperbolicGenerator.h"
#include "../ErdosRenyiGenerator.h"
#include "../ChungLuGenerator.h"
#include "../HavelHakimiGenerator.h"
#include "../RmatGenerator.h"
#include "../../viz/PostscriptWriter.h"
#include "../../community/ClusteringGenerator.h"
#include "../../community/PLP.h"
#include "../../community/PLM.h"
#include "../../io/METISGraphWriter.h"
#include "../../io/DotGraphWriter.h"
#include "../BarabasiAlbertGenerator.h"
#include "../../io/GraphIO.h"
#include "../../io/METISGraphReader.h"
#include "../../properties/GraphProperties.h"
#include "../../community/Modularity.h"
#include "../../dynamics/GraphUpdater.h"
#include "../../auxiliary/MissingMath.h"

namespace NetworKit {

class GeneratorsGTest: public testing::Test {
public:
	GeneratorsGTest();
	virtual ~GeneratorsGTest();

	vector<double> getAngles(DynamicHyperbolicGenerator dynGen) {
		return dynGen.angles;
	}

	vector<double> getRadii(DynamicHyperbolicGenerator dynGen) {
		return dynGen.radii;
	}

	count getQuadTreeHeight(DynamicHyperbolicGenerator dynGen) {
		return dynGen.quadTreeHeight();
	}
};

} /* namespace NetworKit */
#endif /* GENERATORSGTEST_H_ */

#endif /*NOGTEST */
