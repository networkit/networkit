/*
 * LAMGGTest.cpp
 *
 *  Created on: 20.11.2014
 *      Author: Michael
 */

#include "LAMGGTest.h"
#include "../LAMG/MultiLevelSetup.h"
#include "../LAMG/SolverLamg.h"
#include "../../io/LineFileReader.h"
#include "../../auxiliary/Timer.h"

#include "../NearlyLinearLaplacianSmoother.h"
#include "../GaussSeidelRelaxation.h"

namespace NetworKit {

TEST_F(LAMGGTest, testSmallGraphs) {
	METISGraphReader reader;
	GaussSeidelRelaxation gaussSmoother;
	Smoother *smoother = new GaussSeidelRelaxation();
	MultiLevelSetup setup(gaussSmoother);
	Aux::Timer timer;
	for (index i = 0; i < GRAPH_INSTANCES.size(); ++i) {
		string graph = GRAPH_INSTANCES[i];
		Graph G = reader.read("instances/" + graph);
		ConnectedComponents con(G);
		con.run();
		Partition comps = con.getPartition();
		if (comps.numberOfSubsets() > 1) { // disconnected graphs are currently not supported
			continue;
		}

		LevelHierarchy hierarchy;
		timer.start();
		setup.setup(G, hierarchy);
		SolverLamg solver(hierarchy, *smoother);
		timer.stop();
		INFO("setup time\t ", timer.elapsedMilliseconds());

		Vector b = randZeroSum(G, 1234);


//		Vector b(G.numberOfNodes());
//		LineFileReader reader;
//		std::vector<string> lines = reader.read("instances/facebook100/Auburn71.mat.txt_b");
//		std::string input = lines[0];
//		std::string current = "";
//		int idx = 0;
//		for (int i = 0; i < input.size(); ++i) {
//			if (input[i] != ',') {
//				current += input[i];
//			} else {
//				b[idx++] = std::stod(current);
//				current = "";
//			}
//		}
//		b[b.getDimension() - 1] = std::stod(current);

//		INFO("b ", b.length());




		Vector x = randVector(G.upperNodeIdBound(), -1, 1);


//		Vector x(G.numberOfNodes());
//		lines = reader.read("instances/facebook100/Auburn71.mat.txt_x");
//		input = lines[0];
//		current = "";
//		idx = 0;
//		for (int i = 0; i < input.size(); ++i) {
//			if (input[i] != ',') {
//				current += input[i];
//			} else {
//				x[idx++] = std::stod(current);
//				current = "";
//			}
//		}
//		x[x.getDimension() - 1] = std::stod(current);
//
//		INFO("x ", x.length());

//		Vector result(G.upperNodeIdBound());

		LAMGSolverStatus status;
		status.maxConvergenceTime = 10 * 60 * 1000;
		status.desiredResidual = 1e-6;

		Vector result = x;
		INFO("Solving equation system - Gauss-Seidel");
		timer.start();
		solver.solve(result, b, status);
		timer.stop();
		INFO("solve time\t ", timer.elapsedMilliseconds());
		INFO("DONE");

	}

	delete smoother;
}



//TEST_F(LAMGGTest, generateGraphs) {
//	for (count size : grid2DSizes) {
//		genGrid2D(size);
//	}
//
//	for (count size : grid3DSizes) {
//		genGrid3D(size);
//	}
//
//	for (count size : barabasiSizes) {
//		genBarabasi(size);
//	}
//}

Vector LAMGGTest::randVector(count dimension, double lower, double upper) const {
	Vector randVector(dimension);
	for (index i = 0; i < dimension; ++i) {
		randVector[i] = 2.0 * Aux::Random::probability() - 1.0;
	}

	// introduce bias
	for (index i = 0; i < dimension; ++i) {
		randVector[i] = randVector[i] * randVector[i];
	}

	return randVector;
}

Vector LAMGGTest::capacitanceProblem(const Graph &graph) const {
	Vector b(graph.upperNodeIdBound(), 0.0);
	ConnectedComponents con(graph);
	count n = graph.upperNodeIdBound();
	con.run();
	Partition comps = con.getPartition();

	for (index i : comps.getSubsetIds()) {
		std::set<index> members = comps.getMembers(i);
		if (members.size() > 2) {
			count i = 0;
			for (index element : members) {
				if (i == 0) {
					b[element] = 1;
					i++;
				}

				if (i == 1) {
					b[element] = -1;
					break;
				}
			}

			break;
		}
	}

	return b;
}

Vector LAMGGTest::randZeroSum(const Graph& G, size_t seed) const {
	mt19937 rand(seed);
	auto rand_value = uniform_real_distribution<double>(-1.0, 1.0);
	ConnectedComponents con(G);
	count n = G.numberOfNodes();
	con.run();
	Partition comps = con.getPartition();

	/* Fill each component randomly such that its sum is 0 */
	Vector b(n, 0.0);

	for (int id : comps.getSubsetIds()) {
		auto indexes = comps.getMembers(id);
		assert(!indexes.empty());
		double sum = 0.0;
		for (auto entry : indexes) {
			b[entry] = rand_value(rand);
			sum += b[entry];
		}
		b[*indexes.begin()] -= sum;
	}

	return b;
}

// 2D-Grid with unit weights
void LAMGGTest::genGrid2D(count n) const {
  Graph G(n*n);
  for (index i = 0; i < n; ++i) {
    for (index j = 0; j < n; ++j) {
      if (i < n-1) {
        G.addEdge(i*n + j, (i+1)*n + j);
      }
      if (j < n-1) {
        G.addEdge(i*n + j, i*n + (j+1));
      }
    }
  }

  METISGraphWriter writer;
  writer.write(G, false, Aux::toStringF("instances/grid/Laplace_%sx%s.graph", n, n));
}

// 3D-Grid with unit weights
void LAMGGTest::genGrid3D(count n) const {
  Graph G(n*n*n);
  for (index i = 0; i < n; ++i) {
    for (index j = 0; j < n; ++j) {
      for (index k = 0; k < n; ++k) {
        if (i < n-1) {
          G.addEdge(i*n*n + j*n + k, (i+1)*n*n + j*n + k);
        }
        if (j < n-1) {
          G.addEdge(i*n*n + j*n + k, i*n*n + (j+1)*n + k);
        }
        if (k < n-1) {
          G.addEdge(i*n*n + j*n + k, i*n*n + j*n + k+1);
        }
      }
    }
  }

  METISGraphWriter writer;
  writer.write(G, false, Aux::toStringF("instances/grid3/Laplace_%sx%sx%s.graph", n, n, n));
}

/* Preferential attachment random graph with random weights */
void LAMGGTest::genBarabasi(count n, count attachment) const {
  random_device rd;
  mt19937 engine(rd());
  auto rand_weight = uniform_real_distribution<edgeweight>(0.1, 10.0);

  BarabasiAlbertGenerator gen(attachment /* degree */, n, attachment);
  Graph G = Graph(gen.generate(), true, false);
  G.forEdges([&] (node u, node v) {
    G.setWeight(u, v, rand_weight(engine));
  });

  METISGraphWriter writer;
  writer.write(G, true, Aux::toStringF("instances/barabasi/%s_att_%s_unweighted.graph", n, attachment));
}

} /* namespace NetworKit */
