/*
 * ApproxBetweennessGTest.cpp
 *
 *  Created on: 30.06.2014
 *      Author: moritzl
 */

#include "ApproxBetweennessGTest.h"
#include "../ApproxBetweenness.h"
#include "../Betweenness.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../generators/DorogovtsevMendesGenerator.h"
#include "../../properties/Diameter.h"

namespace NetworKit {


TEST_F(ApproxBetweennessGTest, benchApproxDiameterErdos) {
	ErdosRenyiGenerator gen(10000,0.001);
	Graph G1 = gen.generate();
	ApproxBetweenness approx(G1, 0.05, 0.1, 20);
	approx.run();
}

}
