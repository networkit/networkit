/*
 * SSSPGTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include "SSSPGTest.h"
#include "../DynBFS.h"
#include "../BFS.h"
#include "../DirOptBFS.h"
#include "../DynDijkstra.h"
#include "../Dijkstra.h"
#include "../../io/METISGraphReader.h"
#include "../../io/KONECTGraphReader.h"
#include "../../auxiliary/Log.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Random.h"
#include "../../auxiliary/Parallelism.h"
#include <stack>
#include <string>

namespace NetworKit {

TEST_F(SSSPGTest, testDijkstra) {
/* Graph:
         ______
		/      \
	   0    3   6
		\  / \ /
		 2    5
		/  \ / \
	   1    4   7
*/
	int n = 8;
	Graph G(n, true);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);
	G.addEdge(0, 6);


	Dijkstra sssp(G, 5, true, true);
	sssp.run();
	std::vector<node> stack = sssp.getStack();
	while (!stack.empty()) {
		node t = stack.back();
		stack.pop_back();
		DEBUG(t);
	}
}

TEST_F(SSSPGTest, testShortestPaths) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");
	INFO("The graph has been read.");
	int source = 2;
	BFS bfs(G, source);
	bfs.run();
	bigfloat max = 0;
	node x;
	G.forNodes([&](node n){
		if(bfs.numberOfPaths(n) > max) {
			max = bfs.numberOfPaths(n);
			x = n;
		}
	});
	count dist = bfs.distance(x);
	std::set<std::vector<node>> paths = bfs.getPaths(x, true);
	count i = 0;
	for (auto path : paths) {
		INFO("Path number ", i);
		i ++;
		INFO(path);
		EXPECT_EQ(path[0], source);
		EXPECT_EQ(path[dist], x);
	}
	INFO("Maximum number of shortest paths: ", bfs.numberOfPaths(x));
	INFO("Distance: ", dist);
}

TEST_F(SSSPGTest, testGetAllShortestPaths) {
/* Graph:

	   0    3   6   9
		\  / \ / \ /
         2    5   8
		/  \ / \ / \
	   1    4   7   10
*/
	int n = 11;
	Graph G(n, false);
	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);
	G.addEdge(6, 8);
	G.addEdge(7, 8);
	G.addEdge(8, 9);
	G.addEdge(8, 10);
	Dijkstra sssp(G, 0, true, false);
	sssp.run();
	std::set<std::vector<node>> paths = sssp.getPaths(9, true);
	count i = 0;
	for (auto path : paths) {
		INFO("Path number ", i);
		i ++;
		for (node n : path) {
			INFO(n);
		}
	}
}

TEST_F(SSSPGTest, testDirectedBFS) {
/* Graph:
         ________
		/        \.
	   0     3.    6
		\. ./  \ ./
		  2     .5
		./  \. / \.
	   1     4    7
*/
	int n = 8;
	// G directed unweighted
	Graph G(n, false, true);

	G.addEdge(0, 6);
	G.addEdge(0, 2);
	G.addEdge(3, 2);
	G.addEdge(5, 3);
	G.addEdge(6, 5);
	G.addEdge(5, 7);
	G.addEdge(4, 5);
	G.addEdge(2, 4);
	G.addEdge(2, 1);


	BFS sssp(G, 0);
	sssp.run();
	EXPECT_EQ(sssp.distance(0), 0);
	EXPECT_EQ(sssp.distance(1), 2);
	EXPECT_EQ(sssp.distance(2), 1);
	EXPECT_EQ(sssp.distance(3), 3);
	EXPECT_EQ(sssp.distance(4), 2);
	EXPECT_EQ(sssp.distance(5), 2);
	EXPECT_EQ(sssp.distance(6), 1);
	EXPECT_EQ(sssp.distance(7), 3);
}

TEST_F(SSSPGTest, testDirectedDijkstra) {
/* Graph:
         ________
		/        \.
	   0     3.    6
		\. ./  \ ./
		  2     .5
		./  \. / \.
	   1     4    7
*/
	int n = 8;
	// G directed unweighted
	Graph G(n, false, true);

	G.addEdge(0, 6, 1);
	G.addEdge(0, 2, 1);
	G.addEdge(3, 2, 1);
	G.addEdge(5, 3, 1);
	G.addEdge(6, 5, 1);
	G.addEdge(5, 7, 1);
	G.addEdge(4, 5, 1);
	G.addEdge(2, 4, 1);
	G.addEdge(2, 1, 1);


	Dijkstra sssp(G, 0);
	sssp.run();
	EXPECT_EQ(sssp.distance(0), 0);
	EXPECT_EQ(sssp.distance(1), 2);
	EXPECT_EQ(sssp.distance(2), 1);
	EXPECT_EQ(sssp.distance(3), 3);
	EXPECT_EQ(sssp.distance(4), 2);
	EXPECT_EQ(sssp.distance(5), 2);
	EXPECT_EQ(sssp.distance(6), 1);
	EXPECT_EQ(sssp.distance(7), 3);
}

TEST_F(SSSPGTest, testDirOptBFS) {
/* Graph:
         ______
		/      \
	   0    3   6
		\  / \ /
		 2    5
		/  \ / \
	   1    4   7
*/
	int n = 8;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);
	G.addEdge(5, 6);
	G.addEdge(5, 7);
	G.addEdge(0, 6);
	G.indexEdges();

	for (unsigned int storePaths = 0; storePaths < 2; ++storePaths) {
		for (unsigned int storeStack = 0; storeStack < 2; ++storeStack) {

			BFS bfs_ref(G, 5, storePaths, storeStack);
			bfs_ref.run();

			DirOptBFS bfs_diropt(G, 5, storePaths, storeStack);
			bfs_diropt.run();

			auto ref_distances = bfs_ref.getDistances();
			auto do_distances = bfs_diropt.getDistances();


			// computed distances from both BFS should be the same
			// tests: SSSP.getDistances();
			EXPECT_EQ(ref_distances, do_distances) << "distances from source are supposed to be the same";

			if (storeStack) {
				auto ref_stack = bfs_ref.getStack();
				auto do_stack = bfs_diropt.getStack();
				EXPECT_EQ(ref_stack.size(),do_stack.size());
				auto &min_stack = (ref_stack.size() < do_stack.size())?ref_stack:do_stack;
				// the exact order of the stack is unlikely to be the same
				// however, the distance of stack.top() should be the same in both BFS
				// tests: SSSP.getStack();
				while (!min_stack.empty()) {
					auto ref_top = ref_stack.back();
					auto do_top = do_stack.back();
					EXPECT_EQ(ref_distances[ref_top],do_distances[do_top]);
					ref_stack.pop_back(); do_stack.pop_back();
				}
			}

			if (storePaths) {
				// for both BFS, the number of paths for each node should be the same
				// tests: SSSP.numberOfPaths(t)
				count counter = 0;
				count counter_smaller = 0;
				G.forNodes([&](node v){
					//EXPECT_EQ(bfs_ref.numberOfPaths(v),bfs_diropt.numberOfPaths(v)) << "number of paths for node " << v << " differ";
					counter += bfs_ref.numberOfPaths(v) != bfs_diropt.numberOfPaths(v);
					counter_smaller += bfs_ref.numberOfPaths(v) > bfs_diropt.numberOfPaths(v);
				});
				EXPECT_EQ(0,counter) << "of " << G.numberOfNodes() << " are wrong";
				EXPECT_EQ(0,counter_smaller) << "are smaller than the actual number of paths";
			}
		}
	}
}

TEST_F(SSSPGTest, testDirOptBFSOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");
	G.indexEdges();

	for (unsigned int storePaths = 0; storePaths < 2; ++storePaths) {
		for (unsigned int storeStack = 0; storeStack < 2; ++storeStack) {
			BFS bfs_ref(G, 5, storePaths, storeStack);
			bfs_ref.run();

			DirOptBFS bfs_diropt(G, 5, storePaths, storeStack);
			bfs_diropt.run();

			auto ref_distances = bfs_ref.getDistances();
			auto do_distances = bfs_diropt.getDistances();


			// computed distances from both BFS should be the same
			// tests: SSSP.getDistances();
			EXPECT_EQ(ref_distances, do_distances) << "distances from source are supposed to be the same";

			if (storeStack) {
				auto ref_stack = bfs_ref.getStack();
				auto do_stack = bfs_diropt.getStack();
				auto &min_stack = (ref_stack.size() < do_stack.size())?ref_stack:do_stack;
				// the exact order of the stack is unlikely to be the same
				// however, the distance of stack.top() should be the same in both BFS
				// tests: SSSP.getStack();
				while (!min_stack.empty()) {
					auto ref_top = ref_stack.back();
					auto do_top = do_stack.back();
					EXPECT_EQ(ref_distances[ref_top],do_distances[do_top]);
					ref_stack.pop_back(); do_stack.pop_back();
				}
			}

			if (storePaths) {
				// for both BFS, the number of paths for each node should be the same
				// tests: SSSP.numberOfPaths(t)
				count counter = 0;
				count counter_smaller = 0;
				G.forNodes([&](node v){
					//EXPECT_EQ(bfs_ref.numberOfPaths(v),bfs_diropt.numberOfPaths(v)) << "number of paths for node " << v << " differ";
					counter += bfs_ref.numberOfPaths(v) != bfs_diropt.numberOfPaths(v);
					counter_smaller += bfs_ref.numberOfPaths(v) > bfs_diropt.numberOfPaths(v);
				});
				EXPECT_EQ(0,counter) << "of " << G.numberOfNodes() << " are wrong";
				EXPECT_EQ(0,counter_smaller) << "are smaller than the actual number of paths";
			}
		}
	}
}

TEST_F(SSSPGTest, testDirOptBFSOnDirectedRealGraph) {
	KONECTGraphReader reader(' ');
	Graph G = reader.read("input/foodweb-baydry.konect");
	G.indexEdges();

	for (unsigned int storePaths = 0; storePaths < 2; ++storePaths) {
		for (unsigned int storeStack = 0; storeStack < 2; ++storeStack) {
			BFS bfs_ref(G, 5, storePaths, storeStack);
			bfs_ref.run();

			DirOptBFS bfs_diropt(G, 5, storePaths, storeStack);
			bfs_diropt.run();

			auto ref_distances = bfs_ref.getDistances();
			auto do_distances = bfs_diropt.getDistances();


			// computed distances from both BFS should be the same
			// tests: SSSP.getDistances();
			EXPECT_EQ(ref_distances, do_distances) << "distances from source are supposed to be the same";

			if (storeStack) {
				auto ref_stack = bfs_ref.getStack();
				auto do_stack = bfs_diropt.getStack();
				auto &min_stack = (ref_stack.size() < do_stack.size())?ref_stack:do_stack;
				// the exact order of the stack is unlikely to be the same
				// however, the distance of stack.top() should be the same in both BFS
				// tests: SSSP.getStack();
				while (!min_stack.empty()) {
					auto ref_top = ref_stack.back();
					auto do_top = do_stack.back();
					EXPECT_EQ(ref_distances[ref_top],do_distances[do_top]);
					ref_stack.pop_back(); do_stack.pop_back();
				}
			}

			if (storePaths) {
				// for both BFS, the number of paths for each node should be the same
				// tests: SSSP.numberOfPaths(t)
				count counter = 0;
				count counter_smaller = 0;
				G.forNodes([&](node v){
					//EXPECT_EQ(bfs_ref.numberOfPaths(v),bfs_diropt.numberOfPaths(v)) << "number of paths for node " << v << " differ";
					counter += bfs_ref.numberOfPaths(v) != bfs_diropt.numberOfPaths(v);
					counter_smaller += bfs_ref.numberOfPaths(v) > bfs_diropt.numberOfPaths(v);
				});
				EXPECT_EQ(0,counter) << "of " << G.numberOfNodes() << " are wrong";
				EXPECT_EQ(0,counter_smaller) << "are smaller than the actual number of paths";
			}
		}
	}
}

TEST_F(SSSPGTest, benchDirOptBFS) {
	auto minimal_bfs = [](const Graph& G, count source) {
		count z = G.upperNodeIdBound();
		std::vector<double> distances(z,std::numeric_limits<edgeweight>::max());
		std::vector<bool> visited(z,false);
		distances[source] = 0;
		visited[source] = true;
		std::queue<node> q;
		q.push(source);
		while(!q.empty()) {
			auto current = q.front();
			q.pop();
			G.forNeighborsOf(current,[&](node v) {
				if (!visited[v]) {
					visited[v] = true;
					q.push(v);
					distances[v] = distances[current] + 1;
				}
			});
		}
		return std::move(distances);
	};

	bool local = false;
	std::string base;
	std::vector<std::string> datasets;
	if (!local) {
		base = "/algoDaten/staudt/Graphs/Collections/NwkBenchmark/";
		datasets = { "caidaRouterLevel.metis.graph", "coAuthorsDBLP.metis.graph", "in-2004.metis.graph", "con-fiber_big.metis.graph", "uk-2002.metis.graph" };
	} else {
		base = "input/";
		datasets = {"PGPgiantcompo.graph","astro-ph.graph", "wing.graph", "caidaRouterLevel.graph","as-Skitter.metis.graph"};
	}
	METISGraphReader reader;
	Aux::Timer t;

	Aux::Random::setSeed(42, false);
	for (auto& file : datasets) {
		Graph G = reader.read(base+file);
		G.indexEdges();
		count nRuns = 10;
		std::cout << "benchmarking BFS variants: " << nRuns << " runs on " << G.toString() << ", reporting average time in ms" << std::endl;
		std::vector<node> startNodes(nRuns);
		// generate startNode sequence
		for (index i = 0; i < nRuns; ++i) {
			startNodes[i] = Aux::Random::integer(0,G.numberOfNodes());
		}
		t.start();
		for (index i = 0; i < nRuns; ++i) {
			std::ignore = minimal_bfs(G,startNodes[i]);
		}
		t.stop();
		count avg_time_minbfs = t.elapsedMilliseconds() / nRuns;
		std::cout << "minimal bfs:\t" << avg_time_minbfs << std::endl;

		for (int k = 0; k < 2; ++k) {
			for (int j = 0; j < 2; ++j) {
				t.start();
				for (index i = 0; i < nRuns; ++i) {
					BFS bfs(G,startNodes[i],k,j);
					bfs.run();
				}
				t.stop();
				count avg_time_refbfs = t.elapsedMilliseconds() / nRuns;
				std::cout << "reference bfs(storePaths=" << k << ",storeStack=" << j <<"):\t" << avg_time_refbfs << std::endl;

				t.start();
				for (index i = 0; i < nRuns; ++i) {
					DirOptBFS bfs(G,startNodes[i],k,j);
					bfs.run();
				}
				t.stop();
				count avg_time_dobfs = t.elapsedMilliseconds() / nRuns;
				std::cout << "diropt bfs(storePaths=" << k << ",storeStack=" << j <<"):\t" << avg_time_dobfs << std::endl;
			}
		}
		std::cout << "----------------------------------------------" << std::endl;
	}
}

TEST_F(SSSPGTest, benchDirOptBFSThreading) {
	std::vector<std::string> datasets = {
		"input/PGPgiantcompo.graph",
		"input/astro-ph.graph",
		"input/caidaRouterLevel.graph",
		"/algoDaten/staudt/Graphs/Collections/NwkBenchmark/in-2004.metis.graph",
		"/algoDaten/staudt/Graphs/Collections/NwkBenchmark/con-fiber_big.metis.graph",
		//"/algoDaten/staudt/Graphs/Collections/NwkBenchmark/uk-2002.metis.graph"
		"/algoDaten/staudt/Graphs/Static/DIMACS/XLarge/uk-2007-05.metis.graph"
	};
	METISGraphReader reader;
	Aux::Timer t;
	Aux::Random::setSeed(42, false);

	for (auto& file : datasets) {
		Graph G = reader.read(file);
		G.indexEdges();
		count nRuns = 10;
		std::cout << "benchmarking BFS variants: " << nRuns << " runs on " << G.toString() << ", reporting average time in ms" << std::endl;
		std::vector<node> startNodes(nRuns);
		// generate startNode sequence
		for (index i = 0; i < nRuns; ++i) {
			startNodes[i] = Aux::Random::integer(0,G.numberOfNodes());
		}

		for (int k = 0; k < 2; ++k) {
			for (int j = 0; j < 2; ++j) {
				t.start();
				for (index i = 0; i < nRuns; ++i) {
					BFS bfs(G,startNodes[i],k,j);
					bfs.run();
				}
				t.stop();
				count avg_time_refbfs = t.elapsedMilliseconds() / nRuns;
				std::cout << "reference bfs(storePaths=" << k << ",storeStack=" << j <<"):\t" << avg_time_refbfs << std::endl;

				auto threads = {1,2,4,8,16,32,64};
				for (auto& current : threads) {
					Aux::setNumberOfThreads(current);
					t.start();
					for (index i = 0; i < nRuns; ++i) {
						DirOptBFS bfs(G,startNodes[i],k,j);
						bfs.run();
					}
					t.stop();
					count avg_time_dobfs = t.elapsedMilliseconds() / nRuns;
					std::cout << "diropt bfs(storePaths=" << k << ",storeStack=" << j << ",threads=" << current << "):\t" << avg_time_dobfs << std::endl;
				}
			}
		}
		std::cout << "----------------------------------------------" << std::endl;
	}
}



} /* namespace NetworKit */
