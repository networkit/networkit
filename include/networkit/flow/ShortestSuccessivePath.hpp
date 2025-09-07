/*  ShortestSuccessivePath.hpp
*
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef SHORTESTSUCCESSIVEPATH_H
#define SHORTESTSUCCESSIVEPATH_H

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class MinFlowShortestSuccessivePath : public Algorithm {
  using cost = edgeweight;
   public:
        /**
         * @param G            underlying graph (must be directed + weighted)
         * @param capacityAttr name of the edge attribute holding capacity values
         * @param supplyAttr   name of the node attribute holding supply values
         *
         * @throws std::runtime_error if G is not directed or not weighted,
         *         or if either of the named attributes is missing.
         */
        MinFlowShortestSuccessivePath(const Graph &G,
                                      const std::string &capacityName,
                                      const std::string &supplyName);


        void run() override;
        double getTotalCost() const {
          assureFinished();
          return totalCost;
        }
    private:
        const Graph *graph;
        std::string capacityAttributeName;
        std::string supplyAttributeName;
        Graph residualGraph;
        static constexpr const char * FLOW = "flow";
        cost totalCost{};
};

}
#endif //SHORTESTSUCCESSIVEPATH_H
