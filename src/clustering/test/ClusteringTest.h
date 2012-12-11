/*
 * ClusteringTest.h
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#ifndef ClusteringTest_H_
#define ClusteringTest_H_

#include <iostream>

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../aux/log.h"
#include "../Clustering.h"
#include "../ClusteringGenerator.h"
#include "../Modularity.h"

namespace EnsembleClustering {

class ClusteringTest : public CPPUNIT_NS::TestFixture {

	CPPUNIT_TEST_SUITE(ClusteringTest);
	CPPUNIT_TEST(testModularity);
	CPPUNIT_TEST_SUITE_END();

protected:

	GraphGenerator gen;
	Graph randomGraph;


public:

	ClusteringTest();

	virtual ~ClusteringTest();

	void setUp();

	void tearDown();

	/*** Tests **/

	void testModularity();
};

} /* namespace EnsembleClustering */
#endif /* ClusteringTest_H_ */
