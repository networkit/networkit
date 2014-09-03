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

ApproxBetweennessGTest::ApproxBetweennessGTest() {
	// TODO Auto-generated constructor stub

}

ApproxBetweennessGTest::~ApproxBetweennessGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(ApproxBetweennessGTest, testApproxDiameterErdos) {
	ErdosRenyiGenerator gen(10000,0.001);
	Graph G1 = gen.generate();
	ApproxBetweenness approx(G1, 0.05, 0.1, 20);
	approx.run();
}

TEST_F(ApproxBetweennessGTest, testApproxBetweennessSamples) {
	int n = 20000;
	DorogovtsevMendesGenerator generator(n);
	Graph G = generator.generate();
	double epsilon1 = 0.1;
	double epsilon2 = 0.2;
	double epsilon3 = 0.4;
	ApproxBetweenness approx1(G, epsilon1, 0.1);
	ApproxBetweenness approx2(G, epsilon2, 0.1);
	ApproxBetweenness approx3(G, epsilon3, 0.1);
	Betweenness bc(G);
	bc.run();
	approx1.run();
	approx2.run();
	approx3.run();
	std::vector<double> bc_scores = bc.scores();
	std::vector<double> approx_scores1 = approx1.scores();
	std::vector<double> approx_scores2 = approx2.scores();
	std::vector<double> approx_scores3 = approx3.scores();

	double m1 = 0;
	double m2 = 0;
	double m3 = 0;
	G.forNodes([&](node u){
		double norm_value = bc_scores[u]/(n*(n-1));
		if (norm_value-approx_scores1[u] > m1 || norm_value-approx_scores1[u] < -1*m1) {
			m1 = norm_value - approx_scores1[u];
		}
		if (norm_value-approx_scores2[u] > m2 || norm_value-approx_scores2[u] < -1*m2) {
			m2 = norm_value - approx_scores2[u];
		}
		if (norm_value-approx_scores1[u] > m3 || norm_value-approx_scores3[u] < -1*m3) {
			m3 = norm_value - approx_scores3[u];
		}
	});
	count inv1 = 0;
	count inv2 = 0;
	count inv3 = 0;
	G.forNodes([&](node u){
		G.forNodes([&](node v){
			if ((bc_scores[u]-bc_scores[v])*(approx_scores1[u]-approx_scores1[v]) < 0 || (bc_scores[u]-bc_scores[v]!=0 && approx_scores1[u]-approx_scores1[v] == 0))
				inv1++;
			if ((bc_scores[u]-bc_scores[v])*(approx_scores2[u]-approx_scores2[v]) < 0 || (bc_scores[u]-bc_scores[v]!=0 && approx_scores2[u]-approx_scores2[v] == 0))
				inv2++;
			if ((bc_scores[u]-bc_scores[v])*(approx_scores3[u]-approx_scores3[v]) < 0 || (bc_scores[u]-bc_scores[v]!=0 && approx_scores3[u]-approx_scores3[v] == 0))
				inv3++;
		});
	});
	INFO("EPSILON = ", epsilon1);
	INFO("Maximum error: ",m1);
	INFO("Number of inversions: ", inv1);
	INFO("Percentage of inversions: ", double(inv1)/(n*n));
	INFO("Number of samples: ", approx1.numberOfSamples());
	INFO("EPSILON = ", epsilon2);
	INFO("Maximum error: ",m2);
	INFO("Number of inversions: ", inv2);
	INFO("Percentage of inversions: ", double(inv2)/(n*n));
	INFO("Number of samples: ", approx2.numberOfSamples());
	INFO("EPSILON = ", epsilon3);
	INFO("Maximum error: ",m3);
	INFO("Number of inversions: ", inv3);
	INFO("Percentage of inversions: ", double(inv3)/(n*n));
	INFO("Number of samples: ", approx3.numberOfSamples());
}

}
