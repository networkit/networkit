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

   public:
        /**
         * @param graph            underlying graph (must be directed + weighted)
         * @param capacityAttr name of the edge attribute holding capacity values
         * @param supplyAttr   name of the node attribute holding supply values
         *
         * @throws std::runtime_error if G is not directed or not weighted,
         *         or if either of the named attributes is missing.
         */
        MinFlowShortestSuccessivePath(const Graph &graph,
                                      const std::string &capacityName,
                                      const std::string &supplyName);


        void run() override;
    private:
        void initializeResidualGraph();
        std::string capacityAttributeName;
        std::string supplyAttributeName;
        Graph residualGraph;
};

}
#endif //SHORTESTSUCCESSIVEPATH_H
