/*
 * BarabasiAlbertGenerator.h
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#ifndef BarabasiAlbertGenerator_H_
#define BarabasiAlbertGenerator_H_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Generates a scale-free graph using the Barabasi-Albert preferential attachment model.
 */
class BarabasiAlbertGenerator: public StaticGraphGenerator {
private:
    Graph initGraph;
    count k{0}; //!< Attachments made per node
    count nMax{0}; //!< The maximal number of nodes attached
    count n0{0}; //!< The number of initial connected nodes
    bool batagelj; //!< Specifies whether to use batagelj's method or the original one

    /**
     * Implementation of ALG 5 of Batagelj, Brandes: Efficient Generation of Large Random Networks
     * https://kops.uni-konstanz.de/bitstream/handle/123456789/5799/random.pdf?sequence=1
     * Running time is O(n+m)
     * @return The generated graph
     */
    Graph generateBatagelj();

    Graph generateOriginal();

public:
    BarabasiAlbertGenerator() = default;

    /**
     * This generator implements the preferential attachment model as introduced by Barabasi and Albert[1].
     * The original algorithm is very slow and thus, the much faster method from Batagelj and Brandes[2] is
     * implemented and the current default.
     * The original method can be chosen by setting \p batagelj to false.
     * [1] Barabasi, Albert: Emergence of Scaling in Random Networks
     *   http://arxiv.org/pdf/cond-mat/9910332.pdf
     * [2] ALG 5 of Batagelj, Brandes: Efficient Generation of Large Random Networks
     *   https://kops.uni-konstanz.de/bitstream/handle/123456789/5799/random.pdf?sequence=1
     *
     * The generator will emit a simple graph (given that the seed graph is simple), where all
     * new nodes are initially connected to k random neighbors.
     *
     * @param k     Number of attachments per node
     * @param nMax  Number of nodes in the graph
     * @param n0    Number of connected nodes to begin with
     * @param turbo Specifies whether to use Batagelj and Brandes's method rather than the naive one; default: true
     */
    BarabasiAlbertGenerator(count k, count nMax, count n0 = 0, bool turbo=true);

    BarabasiAlbertGenerator(count k, count nMax, const Graph& initGraph, bool turbo=true);

    Graph generate() override;
};

} /* namespace NetworKit */
#endif /* BarabasiAlbertGenerator_H_ */
