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
	LFRGenerator(count n);

	void setDegreeSequence(std::vector<count> degreeSequence);

	void generatePowerlawDegreeSequence(count avgDegree, count maxDegree, double nodeDegreeExp);

	void setCommunitySizeSequence(std::vector<count> communitySizeSequence);

	void setPartition(Partition zeta);

	void generatePowerlawCommunitySizeSequence(count minCommunitySize, count maxCommunitySize, double communitySizeExp);

	void setMu(double mu);

	void setMu(const std::vector<double> & mu);

	void setMuWithBinomialDistribution(double mu);

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

protected:
	virtual std::vector<std::vector<node>> assignNodesToCommunities();
	virtual Graph generateIntraClusterGraph(std::vector< NetworKit::count > intraDegreeSequence, const std::vector< NetworKit::node > &localToGlobalNode);
	virtual Graph generateInterClusterGraph(const std::vector<count> &externalDegreeSequence);
private:
	count n;

	bool hasDegreeSequence;
	std::vector<count> degreeSequence;

	bool hasCommunitySizeSequence;
	std::vector<count> communitySizeSequence;

	bool hasInternalDegreeSequence;
	std::vector<count> internalDegreeSequence;

	bool hasGraph;
	Graph G;

	bool hasPartition;
	Partition zeta;
};

} // namespace NetworKit

#endif // LFRGENERATOR_H
