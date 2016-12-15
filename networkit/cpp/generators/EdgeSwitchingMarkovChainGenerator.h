#ifndef EDGESWITCHINGMARKOVCHAINGENERATOR_H
#define EDGESWITCHINGMARKOVCHAINGENERATOR_H

#include "StaticDegreeSequenceGenerator.h"

namespace NetworKit {

/**
 * Graph generator for generating a random simple graph with exactly the given degree sequence based on the Edge-Switching Markov-Chain method.
 *
 * This implementation is based on the paper
 * "Random generation of large connected simple graphs with prescribed degree distribution" by Fabien Viger and Matthieu Latapy,
 * available at http://www-rp.lip6.fr/~latapy/FV/generation.html, however without preserving connectivity (this could later be added as
 * optional feature).
 *
 * The Havel-Hakami generator is used for the initial graph generation, then the Markov-Chain Monte-Carlo algorithm as described and
 * implemented by Fabien Viger and Matthieu Latapy but without the steps for ensuring connectivity is executed. This should lead to a
 * graph that is drawn uniformly at random from all graphs with the given degree sequence.
 *
 * Note that at most 10 times the number of edges edge swaps are performed (same number as in the abovementioned implementation) and
 * in order to limit the running time, at most 200 times as many attempts to perform an edge swap are made (as certain degree distributions
 * do not allow edge swaps at all).
 */
class EdgeSwitchingMarkovChainGenerator : public StaticDegreeSequenceGenerator {
public:
	/**
	 * Initializes the configuration model generator with the given degree sequence @a sequence.
	 *
	 * If the sequence cannot be realized, optionally if @a ignoreIfRealizable is true, a graph with a different degree
	 * sequence is generated where certain nodes might not have their full degree.
	 *
	 * @param sequence The wanted degree sequence.
	 * @param ignoreIfRealizable Ignore if the sequence is realizable and generate a graph anyway.
	 */
	EdgeSwitchingMarkovChainGenerator(const std::vector< NetworKit::count > &sequence, bool ignoreIfRealizable = false);

	/**
	 * Generate a graph according to the configuration model.
	 *
	 * Issues a INFO log message if the wanted number of edge swaps cannot be performed because of the limit of attempts (see in the description of the class for details).
	 *
	 * @throws std::runtime_error If the sequence cannot be generated and ignoreIfRealizable is false.
	 * @return The generated graph
	 */
	virtual Graph generate() override;
private:
	bool ignoreIfRealizable;
};

} // namespace NetworKit

#endif // EDGESWITCHINGMARKOVCHAINGENERATOR_H
