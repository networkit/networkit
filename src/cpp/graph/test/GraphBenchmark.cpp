/*
 * GraphBenchmark.cpp
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "GraphBenchmark.h"
#include "../NodeMap.h"
#include "../BFS.h"
#include "../BidirectionalBFS.h"
#include "../Dijkstra.h"
#include "../../auxiliary/Random.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"
#include "../../auxiliary/Random.h"
#include "../../auxiliary/Timer.h"
#include <iostream>


namespace NetworKit {

GraphBenchmark::GraphBenchmark() {
	this->n = 1000;
	INFO("n = " , this->n);
}

GraphBenchmark::~GraphBenchmark() {
	// TODO Auto-generated destructor stub
}


// TASK: benchmark edge insertions standard vs raw

TEST_F(GraphBenchmark, edgeInsertions_noop_seq) {
	int64_t n = this->n;
	Aux::Timer runtime;

	Graph G(n);
	int64_t i = 0;
	runtime.start();
	G.forNodePairs([&](node u, node v) {
		i++;
		// G.insertEdge(u, v);
	});
	runtime.stop();

	TRACE("counted i = " , i);

	INFO("[DONE] edgeInsertions_noop_seq (" , runtime.elapsed().count() , " ms)");

}

TEST_F(GraphBenchmark, edgeInsertions_noop_par) {
	int64_t n = this->n;
	Aux::Timer runtime;

	Graph G(n);
	int64_t i = 0;
	runtime.start();
	G.parallelForNodePairs([&](node u, node v) {
		i++;
		// G.insertEdge(u, v);
	});
	runtime.stop();

	TRACE("counted i = " , i);

	INFO("[DONE] edgeInsertions_noop_par (" , runtime.elapsed().count() , " ms)");

}

TEST_F(GraphBenchmark, edgeInsertions_standard_seq) {
	count n = this->n;
	Aux::Timer runtime;

	Graph G(n);
	runtime.start();
	G.forNodePairs([&](node u, node v) {
		G.addEdge(u, v);
	});
	runtime.stop();

	INFO("[DONE] edgeInsertions_standard_seq (" , runtime.elapsed().count() , " ms)");
	EXPECT_EQ((n * (n-1)) / 2, G.numberOfEdges());


}

//TEST_F(GraphBenchmark, edgeInsertions_standard_par) {
//	int64_t n = this->n;
//	Aux::Timer runtime;
//
//	Graph G(n);
//	runtime.start();
//	G.parallelForNodePairs([&](node u, node v) {
//		G.insertEdge(u, v);
//	});
//	runtime.stop();
//
//	INFO("[DONE] edgeInsertions_standard_par(" , runtime.elapsed().count() , " ms)");
//	EXPECT_EQ((n * (n-1)) / 2, G.numberOfEdges());
//
//}
//
//TEST_F(GraphBenchmark, edgeInsertions_raw_seq) {
//	int64_t n = this->n;
//	Aux::Timer runtime;
//
//	Graph G(n);
//	stinger* S = G.asSTINGER();
//
//	runtime.start();
//	for (node u = 1; u <= n; ++u) {
//		for (node v = u + 1; v <= n; ++v) {
//			stinger_insert_edge_pair(S, G.defaultEdgeType, u, v, G.defaultEdgeWeight, G.defaultTimeStamp);
//		}
//	}
//	runtime.stop();
//
//
//	INFO("[DONE] edgeInsertions_raw_seq (" , runtime.elapsed().count() , " ms)");
//	EXPECT_EQ((n * (n-1)) / 2, G.numberOfEdges());
//
//
//}

//TEST_F(GraphBenchmark, edgeInsertions_raw_par) {
//	int64_t n = this->n;
//	Aux::Timer runtime;
//
//	Graph G(n);
//	stinger* S = G.asSTINGER();
//
//	runtime.start();
//	#pragma omp parallel
//	for (node u = 1; u <= n; ++u) {
//		for (node v = u + 1; v <= n; ++v) {
//			stinger_insert_edge_pair(S, G.defaultEdgeType, u, v, G.defaultEdgeWeight, G.defaultTimeStamp);
//		}
//	}
//	runtime.stop();
//
//	INFO("[DONE] edgeInsertions_raw_par (" , runtime.elapsed().count() , " ms)");
//	EXPECT_EQ((n * (n-1)) / 2, G.numberOfEdges());
//
//}




// Task: precompute incident weights with different methods



TEST_F(GraphBenchmark, weightedDegree_standard_seq) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	NodeMap<double> weightedDegree(n, 0.0);

	G.forNodes([&](node v) {
		weightedDegree[v] = G.weightedDegree(v);
	});
	runtime.stop();

	INFO("[DONE] (" , runtime.elapsed().count() , " ms)");

	// test correctness of result
	bool correct = true;
	G.forNodes([&](node v){
		correct &= (weightedDegree[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}


// TEST: use different containers
// RESULT: NodeMap, vector and array are about equally fast


// TEST: parallelize

TEST_F(GraphBenchmark, weightedDegree_standard_par) {
	int64_t n = this->n;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	Aux::Timer runtime;

	runtime.start();
	NodeMap<double> weightedDegree(n, 0.0);

	G.parallelForNodes([&](node v) {
		weightedDegree[v] = G.weightedDegree(v);
	});
	runtime.stop();

	INFO("[DONE] (" , runtime.elapsed().count() , " ms)");

	// test correctness of result
	bool correct = true;
	G.forNodes([&](node v){
		correct &= (weightedDegree[v] == (n - 1));
	});

	EXPECT_TRUE(correct);
}


// RESULT: significant super-linear speedup regardless of target container

//TEST_F(GraphBenchmark, weightedDegree_raw_seq) {
//	int64_t n = this->n;
//	GraphGenerator graphGen;
//	Graph G = graphGen.makeCompleteGraph(n);
//	stinger* S = G.asSTINGER();
//
//	Aux::Timer runtime;
//
//	runtime.start();
//	NodeMap<double> weightedDegree(n, 0.0);
//
//	for (node v = 1; v <= n; ++v) {
//		double iw = 0.0;
//		STINGER_READ_ONLY_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
//			iw += stinger_edgeweight(S, STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, G.defaultEdgeType);
//		} STINGER_READ_ONLY_FORALL_EDGES_OF_VTX_END();
//		weightedDegree[v] = iw;
//	}
//	runtime.stop();
//
//	INFO("[DONE] (" , runtime.elapsed().count() , " ms)");
//
//	// test correctness of result
//	bool correct = true;
//	G.forNodes([&](node v){
//		correct &= (weightedDegree[v] == (n - 1));
//	});
//
//	EXPECT_TRUE(correct);
//
//}

//
//TEST_F(GraphBenchmark, weightedDegree_raw_par) {
//	int64_t n = this->n;
//	GraphGenerator graphGen;
//	Graph G = graphGen.makeCompleteGraph(n);
//	stinger* S = G.asSTINGER();
//
//	Aux::Timer runtime;
//
//	runtime.start();
//	NodeMap<double> weightedDegree(n, 0.0);
//
//	#pragma omp parallel for
//	for (node v = 1; v <= n; ++v) {
//		double iw = 0.0;
//		STINGER_READ_ONLY_FORALL_EDGES_OF_VTX_BEGIN(S, v) {
//			iw += stinger_edgeweight(S, STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, G.defaultEdgeType);
//		} STINGER_READ_ONLY_FORALL_EDGES_OF_VTX_END();
//		weightedDegree[v] = iw;
//	}
//	runtime.stop();
//
//	INFO("[DONE] (" , runtime.elapsed().count() , " ms)");
//
//	// test correctness of result
//	bool correct = true;
//	G.forNodes([&](node v){
//		correct &= (weightedDegree[v] == (n - 1));
//	});
//
//	EXPECT_TRUE(correct);
//
//}


TEST_F(GraphBenchmark, benchBFS) {
	std::vector<std::string> graphs = {"input/PGPgiantcompo.graph",
		"input/astro-ph.graph",
		"input/caidaRouterLevel.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/Large/europe-osm.metis.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/Street/germany.osm.metis.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/Street/asia.osm.metis.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/XLarge/uk-2007-05.metis.graph"
	};

	for(auto file : graphs) {
//		std::ifstream f(file);
		// bench the graph, if the file exists.
//		if (f.good()) {i
		try {
			METISGraphReader reader;
			Graph G = reader.read(file);
			BFS bfs(G,0);
			Aux::Timer timer;
			count runs = 100;
			std::vector<node> targets;
			for (index i = 0; i < runs; ++i) {
				targets.push_back(Aux::Random::integer(0, G.upperNodeIdBound()));
			}
			std::vector<std::vector<node>> paths;
			timer.start();
			for (index i = 0; i < runs; ++i) {
				bfs.run();
				paths.push_back(bfs.getPath(targets[i]));
			}
			timer.stop();
			//INFO(runs, " with BFS.run() on ", G.getName(), " took:\t", timer.elapsedMilliseconds()/100, " ms");
			std::cout << "BFS.run() on " <<  G.getName() << " took:\t" << timer.elapsedMilliseconds()/100 << " ms" << std::endl;

			std::vector<std::vector<node>> pathsUntil;
			timer.start();
			for (index i = 0; i < runs; ++i) {
				bfs.runUntil(targets[i]);
				pathsUntil.push_back(bfs.getPath(targets[i]));
			}
			timer.stop();
//			INFO(runs, " with BFS.runUntil() on ", G.getName(), " took:\t", timer.elapsedMilliseconds()/100, " ms");
			std::cout <<  "BFS.runUntil("<< ") on "<< G.getName()<< " took:\t"<< timer.elapsedMilliseconds()/100 << " ms" <<std::endl;

			for (index i = 0; i < runs; ++i) {
				EXPECT_EQ(paths[i].size(),pathsUntil[i].size());
			}

		} catch (std::exception e) {
			ERROR("couldnt bench on graph ",file," because of:\t",e.what());
		}
//		}
	}
}

TEST_F(GraphBenchmark, benchBidirBFS) {
	std::vector<std::string> graphs = {"input/PGPgiantcompo.graph",
		"input/astro-ph.graph",
		"input/caidaRouterLevel.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/Large/europe-osm.metis.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/Street/germany.osm.metis.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/Street/asia.osm.metis.graph",
//		"/algoDaten/staudt/Graphs/Static/DIMACS/XLarge/uk-2007-05.metis.graph"
	};

	for(auto file : graphs) {
//		std::ifstream f(file);
		// bench the graph, if the file exists.
//		if (f.good()) {i
		try {
			METISGraphReader reader;
			Graph G = reader.read(file);
			BFS bfs(G,0);
			Aux::Timer timer;
			count runs = 100;
			std::vector<node> targets;
			for (index i = 0; i < runs; ++i) {
				targets.push_back(Aux::Random::integer(0, G.upperNodeIdBound()));
			}
			std::vector<std::vector<node>> paths;
			timer.start();
			for (index i = 0; i < runs; ++i) {
				bfs.run();
				paths.push_back(bfs.getPath(targets[i]));
			}
			timer.stop();
			//INFO(runs, " with BFS.run() on ", G.getName(), " took:\t", timer.elapsedMilliseconds()/100, " ms");
			std::cout << "BFS.run() on " <<  G.getName() << " took:\t" << timer.elapsedMilliseconds()/100 << " ms" << std::endl;
			BidirectionalBFS bbfs(G);
			std::vector<std::vector<node>> pathsUntil;
			timer.start();
			for (index i = 0; i < runs; ++i) {
				bbfs.run(0,targets[i]);
				pathsUntil.push_back(bfs.getPath(targets[i]));
			}
			timer.stop();
//			INFO(runs, " with BFS.runUntil() on ", G.getName(), " took:\t", timer.elapsedMilliseconds()/100, " ms");
			std::cout <<  "BFS.runUntil("<< ") on "<< G.getName()<< " took:\t"<< timer.elapsedMilliseconds()/100 << " ms" <<std::endl;

			for (index i = 0; i < runs; ++i) {
				EXPECT_EQ(paths[i].size(),pathsUntil[i].size());
			}

		} catch (std::exception e) {
			ERROR("couldnt bench on graph ",file," because of:\t",e.what());
		}
//		}
	}
}


TEST_F(GraphBenchmark, benchDijkstra) {
	std::vector<std::string> graphs = {"input/PGPgiantcompo.graph",
		"input/astro-ph.graph",
		"input/caidaRouterLevel.graph",
		"/algoDaten/staudt/Graphs/Static/DIMACS/Large/europe-osm.metis.graph",
		"/algoDaten/staudt/Graphs/Static/DIMACS/Street/germany.osm.metis.graph",
		"/algoDaten/staudt/Graphs/Static/DIMACS/Street/asia.osm.metis.graph",
		"/algoDaten/staudt/Graphs/Static/DIMACS/XLarge/uk-2007-05.metis.graph"
	};

	for(auto file : graphs) {
//		std::ifstream f(file);
		// bench the graph, if the file exists.
//		if (f.good()) {i
		try {
			METISGraphReader reader;
			Graph G = reader.read(file);
			if (!G.isWeighted()) {
				WARN("graph ",file," is not weighted, continue with next graph");
				//continue;
			}
			Dijkstra dij(G,0);
			Aux::Timer timer;
			count runs = 100;
			std::vector<node> targets;
			for (index i = 0; i < runs; ++i) {
				targets.push_back(Aux::Random::integer(0, G.upperNodeIdBound()));
			}
			std::vector<std::vector<node>> paths;
			timer.start();
			for (index i = 0; i < runs; ++i) {
				dij.run();
				paths.push_back(dij.getPath(targets[i]));
			}
			timer.stop();
			//INFO(runs, " with BFS.run() on ", G.getName(), " took:\t", timer.elapsedMilliseconds()/100, " ms");
			std::cout << "Dijkstra.run() on " <<  G.getName() << " took:\t" << timer.elapsedMilliseconds()/100 << " ms" << std::endl;

			std::vector<std::vector<node>> pathsUntil;
			timer.start();
			for (index i = 0; i < runs; ++i) {
				dij.runUntil(targets[i]);
				pathsUntil.push_back(dij.getPath(targets[i]));
			}
			timer.stop();
//			INFO(runs, " with BFS.runUntil() on ", G.getName(), " took:\t", timer.elapsedMilliseconds()/100, " ms");
			std::cout <<  "Dijkstra.runUntil("<< ") on "<< G.getName()<< " took:\t"<< timer.elapsedMilliseconds()/100 << " ms" <<std::endl;

			for (index i = 0; i < runs; ++i) {
				EXPECT_EQ(paths[i].size(),pathsUntil[i].size());
			}

		} catch (std::exception e) {
			ERROR("couldnt bench on graph ",file," because of:\t",e.what());
		}
//		}
	}


}


} /* namespace NetworKit */

#endif /*NOGTEST */
