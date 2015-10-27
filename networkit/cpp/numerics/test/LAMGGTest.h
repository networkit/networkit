/*
 * LAMGGTest.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael
 */

#ifndef LAMGGTEST_H_
#define LAMGGTEST_H_

#include "gtest/gtest.h"

#include "../../algebraic/Matrix.h"
#include "../../algebraic/Vector.h"
#include "../../io/METISGraphReader.h"
#include "../../io/METISGraphWriter.h"
#include "../../generators/BarabasiAlbertGenerator.h"
#include "../../properties/ConnectedComponents.h"
#include "../../structures/Partition.h"

using namespace std;

namespace NetworKit {

class LAMGGTest : public testing::Test {
protected:
	const vector<string> GRAPH_INSTANCES = {/*"undirected_graphs/airfoil1.graph", "undirected_graphs/m14b.graph", "2D_3D_Problems/3elt.graph", "2D_3D_Problems/brack2.graph", "2D_3D_Problems/crack.graph", "2D_3D_Problems/wave.graph", "2D_3D_Problems/whitaker3.graph",
			"acoustics_problem/vibrobox.graph", "circuit_simulation_problem/add20.graph", "circuit_simulation_problem/add32.graph", "circuit_simulation_problem/memplus.graph",
			"economic_problem/finan512.graph", "structural_problem/bcsstk29.graph", "structural_problem/bcsstk30.graph", "structural_problem/bcsstk31.graph", "structural_problem/bcsstk32.graph",
			"structural_problem/bcsstk33.graph", "undirected_graphs/598a.graph", "undirected_graphs/cs4.graph", "undirected_graphs/cti.graph", "undirected_graphs/data.graph", "undirected_graphs/fe_4elt2.graph",
			"undirected_graphs/fe_body.graph", "undirected_graphs/fe_ocean.graph", "undirected_graphs/fe_rotor.graph", "undirected_graphs/fe_sphere.graph", "undirected_graphs/fe_tooth.graph",
			"undirected_graphs/t60k.graph", "undirected_graphs/uk.graph", "undirected_graphs/wing_nodal.graph", */"facebook100/Auburn71.mat.txt.metis"/*"grid/Laplace_512x512.graph"*/};

	const std::vector<count> grid2DSizes = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
	const std::vector<count> grid3DSizes = {2, 4, 8, 16, 32, 64, 128};
	const std::vector<count> barabasiSizes = {25, 100, 500, 1000, 5000, 10000, 50000};

	Vector randZeroSum(const Graph &graph, size_t seed) const;
	Vector randVector(count dimension, double lower, double upper) const;
	Vector capacitanceProblem(const Graph &graph) const;

	void genGrid2D(count n) const;
	void genGrid3D(count n) const;
	void genBarabasi(count n, count attachment = 4) const;

};

} /* namespace NetworKit */

#endif /* LAMGGTEST_H_ */
