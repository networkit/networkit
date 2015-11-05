/*
 * SpanningTreeGTest.cpp
 *
 *  Created on: 20.06.2015
 *      Author: Henning
 */

#include "SpanningTreeGTest.h"
#include "../PseudoRandomSpanningTree.h"
#include "../RandomSpanningTree.h"
#include "../../graph/Graph.h"
#include "../../graph/Sampling.h"
#include "../../graph/BFS.h"
#include "../../io/METISGraphReader.h"
#include <cmath>
#include "omp.h"

namespace NetworKit {

TEST_F(SpanningTreeGTest, testRandomSpanningTree) {
	METISGraphReader reader;
	std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

	for (auto graphname: graphs) {
		std::string filename = "input/" + graphname + ".graph";
		Graph G = reader.read(filename);
		RandomSpanningTree rst(G);
		rst.run();
		Graph T = rst.getTree();

		T.forNodes([&](node u) {
			EXPECT_GE(G.degree(u), 0);
		});

		node r1 = Sampling::randomNode(G);
		node r2 = Sampling::randomNode(G);
		while (r1 == r2) {
			r2 = Sampling::randomNode(G);
		}

		BFS bfs(T, r1, false, false, r2);
		bfs.run();
		EXPECT_LE(bfs.distance(r2), G.numberOfNodes() - 1);
	}
}

TEST_F(SpanningTreeGTest, testPseudoRandomSpanningTree) {
  // TODO: see above
}

TEST_F(SpanningTreeGTest, benchRandomSpanningTree) {
	METISGraphReader reader;
	std::vector<std::string> graphs = {"karate", "PGPgiantcompo", "power", "jazz", "celegans_metabolic", "airfoil1"};
	count reps = 500;

	for (auto graphname: graphs) {
		std::string filename = "input/" + graphname + ".graph";
		Graph G = reader.read(filename);

		Graph Gwr(G.numberOfNodes(), true, false);
		G.forEdges([&](node u, node v) {
			Gwr.addEdge(u, v, 0.0);
		});
		Graph Gwp = Gwr;

		// random sampling
		double rstTime = 0.0;
		RandomSpanningTree rst(G);
		for (index i = 0; i < reps; ++i) {
			double time = omp_get_wtime();
			rst.run();
			rstTime += omp_get_wtime() - time;
			Graph tree = rst.getTree();
			tree.forEdges([&](node u, node v) {
				Gwr.setWeight(u, v, 1 + Gwr.weight(u, v));
			});
		}

		// sampling of pseudo random trees
		double prstTime = 0.0;
		PseudoRandomSpanningTree prst(G);
		for (index i = 0; i < reps; ++i) {
			double time = omp_get_wtime();
			prst.run();
			prstTime += omp_get_wtime() - time;
			Graph tree = prst.getTree();
			tree.forEdges([&](node u, node v) {
				Gwp.setWeight(u, v, 1 + Gwp.weight(u, v));
			});
		}

		// sampling results
		double maxDev = 0.0;
		double l1Dev = 0.0;
		double l2Dev = 0.0;

		double maxRatio = 0.0;
		double minRatio = 1e40;
		double gmeanRatio = 1.0;

		G.forEdges([&](node u, node v) {
			double dev = (Gwr.weight(u, v) - Gwp.weight(u, v)) / reps;
			l1Dev += fabs(dev);
			l2Dev += dev * dev;
			if (dev > maxDev) {
				maxDev = dev;
			}

			if (std::min(Gwr.weight(u,v), Gwp.weight(u, v)) > 0.0) {
				double ratio = ((double) Gwr.weight(u, v) / (double) Gwp.weight(u, v));
				gmeanRatio *= ratio;
				if (ratio > maxRatio) {
					maxRatio = ratio;
				}
				if (ratio < minRatio) {
					minRatio = ratio;
				}
			}
		});
		l2Dev = sqrt(l2Dev);
		gmeanRatio = sqrt(gmeanRatio);
		INFO(graphname, ", max: ", maxDev, ", l1: ", l1Dev, ", l2: ", l2Dev);
		INFO(graphname, " ==> time ratio: ", (prstTime / rstTime), ", maxRatio: ", maxRatio, ", minRatio: ", minRatio, ", gmeanRatio: ", gmeanRatio);

		// TODO: ggf. besser als externes Programm
		// TODO: random shuffle gemaess Knotengrad, die also weiter nach vorne
		// mit hoeherer Wkt.
		// Behelf: Kantenarray mit Multiplizitaet jeder Kante gemaess Summe der inzidenten Knoten
		// daraus dann zufaellig ziehen, nur beruecksichtigen, was noch nicht gezogen wurde
		// langsamer, aber nur als proof of concept!
	}
}


} /* namespace NetworKit */
