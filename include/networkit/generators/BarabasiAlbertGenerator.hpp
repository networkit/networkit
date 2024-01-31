/*
 * BarabasiAlbertGenerator.hpp
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#ifndef NETWORKIT_GENERATORS_BARABASI_ALBERT_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_BARABASI_ALBERT_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Generates a scale-free graph using the Barabasi-Albert preferential attachment model.
 */
class BarabasiAlbertGenerator final : public StaticGraphGenerator {
    Graph initGraph;
    count k{0};      //!< Attachments made per node
    count nMax{0};   //!< The maximal number of nodes attached
    count n0{0};     //!< The number of initial connected nodes
    bool sequential; //!< Specifies whether to compute sequentially using batagelj's method or in
                     //!< parallel (if the number of threads allow for it).

    Graph generateParallel();
    Graph generateBatagelj();

public:
    BarabasiAlbertGenerator() = default;

    /**
     * This generator uses the preferential attachment model as introduced by Barabasi and
     * Albert[1], implemented in the much faster method from Batagelj and Brandes[2] per default
     * where the running time is O(n+m). Furthermore there is a parallel version from Sanders and
     * Schulz[3] implemented. This implementation can be selected by setting sequential=false.
     * Empirically, the parallel version shows better runtimes (when executed multi-threaded).
     *
     * [1] Barabasi, Albert: [Emergence of Scaling in Random Networks]
     *   (http://arxiv.org/pdf/cond-mat/9910332.pdf)
     * [2] ALG 5 of Batagelj, Brandes: [Efficient Generation of Large Random Networks]
     *   (https://kops.uni-konstanz.de/bitstream/handle/123456789/5799/random.pdf?sequence=1)
     * [3] Peter Sanders, Christian Schulz: [Scalable generation of scale-free graphs]
     *   (https://www.sciencedirect.com/science/article/pii/S0020019016300102)
     *
     * The generator will emit a simple graph, where all
     * new nodes are initially connected to k random neighbors.
     *
     * @param k     Number of attachments per node
     * @param nMax  Number of nodes in the graph
     * @param n0    Number of connected nodes to begin with
     * @param sequential Specifies whether to use Batagelj and Brandes's method (sequential) or a
     * parallel variant. Default: true.
     */
    BarabasiAlbertGenerator(count k, count nMax, count n0 = 0, bool sequential = true);

    /**
     * This generator uses the preferential attachment model as introduced by Barabasi and
     * Albert[1], implemented in the much faster method from Batagelj and Brandes[2] per default
     * where the running time is O(n+m). Furthermore there is a parallel version from Sanders and
     * Schulz[3] implemented. This implementation can be selected by setting sequential=false.
     * Empirically, the parallel version shows better runtimes (when executed multi-threaded).
     *
     * [1] Barabasi, Albert: [Emergence of Scaling in Random Networks]
     *   (http://arxiv.org/pdf/cond-mat/9910332.pdf)
     * [2] ALG 5 of Batagelj, Brandes: [Efficient Generation of Large Random Networks]
     *   (https://kops.uni-konstanz.de/bitstream/handle/123456789/5799/random.pdf?sequence=1)
     * [3] Peter Sanders, Christian Schulz: [Scalable generation of scale-free graphs]
     *   (https://www.sciencedirect.com/science/article/pii/S0020019016300102)
     *
     * The generator will emit a simple graph, where all
     * new nodes are initially connected to k random neighbors.
     *
     * @param k     Number of attachments per node
     * @param nMax  Number of nodes in the graph
     * @param initGraph   The initial graph to start from
     * @param sequential Specifies whether to use Batagelj and Brandes's method (sequential) or a
     * parallel variant. Default: true.
     */
    BarabasiAlbertGenerator(count k, count nMax, const Graph &initGraph, bool sequential = true);

    Graph generate() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_BARABASI_ALBERT_GENERATOR_HPP_
