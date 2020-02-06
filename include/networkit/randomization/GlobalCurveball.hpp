/*
 * GlobalCurveball.hpp
 *
 *  Created on: 26.05.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
// networkit-format

#ifndef NETWORKIT_RANDOMIZATION_GLOBAL_CURVEBALL_HPP_
#define NETWORKIT_RANDOMIZATION_GLOBAL_CURVEBALL_HPP_

#include <memory>
#include <utility>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// pImpl
namespace CurveballDetails {
class GlobalCurveballImpl;
}

class GlobalCurveball final : public Algorithm {
public:
    /**
     * Instantiate a GlobalCurveball object
     *
     * @param G                        Undirected and unweighted graph to be randomized
     * @param number_of_global_trades  Number of global trades to be executed (each edge
     *                                 is considered exactly twice per global traded)
     * @param allowSelfLoops           May only be set if graph is directed.
     * @param degreePreservingShufflePreprocessing Execute DegreePreservingShuffle
     *                                 (see Algorithm for description) as a preprocessing
     *                                 step. This is more efficient than calling the algorithm
     *                                 explicitly.
     *
     * @note Self loops can only be realized for directed graphs.
     * @warning If self loops are forbidden, degreePreservingShuffle is necessary for
     * directed graphs, since otherwise some topologies cannot be realized (i.e., only
     * preprocessing allows for uniform samples).
     */
    explicit GlobalCurveball(const Graph &G, count number_of_global_trades = 20,
                             bool allowSelfLoops = false,
                             bool degreePreservingShufflePreprocessing = true);

    ~GlobalCurveball();

    /**
     * Execute trades as configured in the constructor.
     * @warning This function has to be called exactly one before invoking getGraph()
     */
    void run() override final;

    /**
     * Returns a new graph instance with the same degree sequence as the input
     * graph, but with randomized neighbourhoods.
     */
    Graph getGraph();

    std::string toString() const override final;

    bool isParallel() const override final { return false; }

private:
    std::unique_ptr<CurveballDetails::GlobalCurveballImpl> impl;
    unsigned numGlobalTrades;
    bool degreePreservingShuffle;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_GLOBAL_CURVEBALL_HPP_
