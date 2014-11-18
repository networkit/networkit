#ifndef LFRGENERATOR_H
#define LFRGENERATOR_H

#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * The LFR clustered graph generator as introduced by Andrea Lancichinetti, Santo Fortunato, and Filippo Radicchi.
 *
 * The community assignment follows the algorithm described in
 * "Benchmark graphs for testing community detection algorithms". The edge generation is however taken from their follow-up publication
 * "Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities". Parts of the
 * implementation follow the choices made in their implementation which is available at https://sites.google.com/site/andrealancichinetti/software
 * but other parts differ, for example some more checks for the realizability of the community and degree size distributions are done
 * instead of heavily modifying the distributions.
 *
 * The configuration model implementation in NetworKit is used which is different from the implementation in the original LFR benchmark.
 */
class LFRGenerator {
public:
	LFRGenerator(count n, count avgDegree, count maxDegree, double mu, double nodeDegreeExp, double communitySizeExp, count minCommunitySize, count maxCommunitySize);

	/**
	 * Generates the graph and the community structure.
	 */
	void run();

	/**
	 * Returns (a copy of) the generated graph.
	 *
	 * @return The generated graph.
	 */
	Graph getGraph() const;

	/**
	 * Returns the generated graph using move semantics.
	 *
	 * @return The generated graph.
	 */
	Graph&& getMoveGraph();

	/**
	 * Returns (a copy of) the generated partition
	 *
	 * @return The generated graph.
	 */
	Partition getPartition() const;

	/**
	 * Returns the generated partition using move semantics.
	 *
	 * @return The generated partition.
	 */
	Partition&& getMovePartition();
private:
	count n, avgDegree, maxDegree;
	double mu, nodeDegreeExp, communitySizeExp;
	count minCommunitySize, maxCommunitySize;

	Graph G;
	Partition zeta;
	bool hasGraph;
	bool hasPartition;
};

} // namespace NetworKit

#endif // LFRGENERATOR_H
