/*
 * DistanceGTest.cpp
 *
 *  Created on: Sep 04, 2015
 *      Author: Maximilian Vogel
 */
#ifndef NOGTEST

#include "DistanceGTest.h"

#include "../Diameter.h"
#include "../EffectiveDiameter.h"

#include "../../generators/DorogovtsevMendesGenerator.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

TEST_F(DistanceGTest, testVertexDiameterPedantically) {
    DorogovtsevMendesGenerator generator(1000);
    Graph G1 = generator.generate();
    Graph G = Graph(G1, true, false);
    count vd = Diameter::estimatedVertexDiameterPedantic(G);
    EXPECT_EQ(1000, vd);
}

TEST_F(DistanceGTest, testExactDiameter) {

   using namespace std;

   vector<pair<string, count>> testInstances= {pair<string, count>("lesmis", 14),
                                               pair<string, count>("jazz", 6),
                                               pair<string, count>("celegans_metabolic", 7)
                                              };

   for (auto testInstance : testInstances) {
       METISGraphReader reader;
       Graph G = reader.read("input/" + testInstance.first + ".graph");
       count diameter = Diameter::exactDiameter(G);
       EXPECT_EQ(diameter, testInstance.second);
   }
}


TEST_F(DistanceGTest, testEstimatedDiameterRange) {

   using namespace std;

   vector<pair<string, count>> testInstances= {
                                               pair<string, count>("celegans_metabolic", 7),
                                               pair<string, count>("jazz", 6)
                                              };

   for (auto testInstance : testInstances) {
       METISGraphReader reader;
       Graph G = reader.read("input/" + testInstance.first + ".graph");
       std::pair<count, count> range = Diameter::estimatedDiameterRange(G, 0.1);
       EXPECT_GE(testInstance.second, range.first);
       EXPECT_LE(testInstance.second, range.second);
   }
}
TEST_F(DistanceGTest, testPedanticDiameterErdos) {
	count n = 5000;
	ErdosRenyiGenerator gen(n,0.001);
	Graph G1 = gen.generate();
	count diameter = Diameter::estimatedVertexDiameterPedantic(G1);
	ASSERT_LE(diameter, n);
}



TEST_F(DistanceGTest, testEffectiveDiameter) {

using namespace std;

vector<string> testInstances= {"celegans_metabolic", "jazz", "lesmis"};

for (auto testInstance : testInstances) {
	METISGraphReader reader;
	Graph G = reader.read("input/" + testInstance + ".graph");
	double effective = EffectiveDiameter::effectiveDiameter(G);
	count exact = Diameter::estimatedDiameterRange(G, 0).first;
	EXPECT_LE(effective, exact);
}
}

TEST_F(DistanceGTest, testEffectiveDiameterExact) {

using namespace std;

vector<string> testInstances= {"celegans_metabolic", "jazz", "lesmis"};

for (auto testInstance : testInstances) {
	METISGraphReader reader;
	Graph G = reader.read("input/" + testInstance + ".graph");
	double effective = EffectiveDiameter::effectiveDiameterExact(G);
	count exact = Diameter::estimatedDiameterRange(G, 0).first;
	EXPECT_LE(effective, exact);
}

const double tol = 1e-3;

/* Graph: n=20, threshold: 20*0.9 = 18 nodes
	1--3--5--7---9
	|  |  |  |   |
	2--4--6--8--10
		|     |
		11----12
			|
		13--14--15
			|
		18--16--17--19
				|
				20
Number of steps needed per node: (1-20)
(7+6+6+5+6+5+5+4+6+5+4+4+5+4+5+5+6+6+7+7) / 20 = 5.4
*/
	count n1 = 20;
	Graph G1(n1);

	G1.addEdge(0,1);
	G1.addEdge(0,2);
	G1.addEdge(1,3);
	G1.addEdge(2,3);
	G1.addEdge(2,4);
	G1.addEdge(3,5);
	G1.addEdge(3,10);
	G1.addEdge(4,5);
	G1.addEdge(4,6);
	G1.addEdge(5,7);
	G1.addEdge(6,8);
	G1.addEdge(6,7);
	G1.addEdge(7,9);
	G1.addEdge(7,11);
	G1.addEdge(8,9);
	G1.addEdge(10,11);
	G1.addEdge(11,13);
	G1.addEdge(12,13);
	G1.addEdge(13,14);
	G1.addEdge(13,15);
	G1.addEdge(15,16);
	G1.addEdge(15,17);
	G1.addEdge(16,18);
	G1.addEdge(16,19);

	double effective1 = EffectiveDiameter::effectiveDiameterExact(G1);
	EXPECT_NEAR(5.4, effective1, tol);

	/* Graph: n=21, threshold: 21*0.9 = 18.9 => 19 nodes
				13---------------3
					|               |
				---14--12--|        |
				|   |   |  |        |
	1--21--18--16--15   |  |        |
				|       |  |        |
		20--17------10--8        |
				|       |  |        |
			19       9--7--5--6--4--11
									|
									2
Number of steps needed per node: (1-21)
(8+7+5+6+6+6+5+5+5+5+7+5+4+4+5+5+5+6+6+6+7) / 21 = 5.619047
*/
	count n2 = 21;
	Graph G2(n2);

	G2.addEdge(0,20);
	G2.addEdge(1,3);
	G2.addEdge(2,3);
	G2.addEdge(2,12);
	G2.addEdge(3,5);
	G2.addEdge(3,10);
	G2.addEdge(4,5);
	G2.addEdge(4,6);
	G2.addEdge(6,7);
	G2.addEdge(6,8);
	G2.addEdge(7,9);
	G2.addEdge(7,11);
	G2.addEdge(8,9);
	G2.addEdge(9,11);
	G2.addEdge(9,16);
	G2.addEdge(11,13);
	G2.addEdge(12,13);
	G2.addEdge(13,14);
	G2.addEdge(13,15);
	G2.addEdge(14,15);
	G2.addEdge(15,16);
	G2.addEdge(15,17);
	G2.addEdge(16,18);
	G2.addEdge(16,19);
	G2.addEdge(17,20);

	double effective2 = EffectiveDiameter::effectiveDiameterExact(G2);
	EXPECT_NEAR(5.619047, effective2, tol);
}

TEST_F(DistanceGTest, testHopPlot) {
	using namespace std;

	vector<string> testInstances= {"celegans_metabolic", "power", "lesmis"};

	const double tol = 1e-2;

	for (auto& testInstance : testInstances) {
		METISGraphReader reader;
		Graph G = reader.read("input/" + testInstance + ".graph");
		map<count, double> hopPlot = EffectiveDiameter::hopPlot(G);
		for (count i=1; i < hopPlot.size(); i++) {
			EXPECT_LE(hopPlot[i-1], hopPlot[i]+tol);
		}
	}
}

} /* namespace NetworKit */

#endif /*NOGTEST */
