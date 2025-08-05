/*  ShortestSuccessivePath.hpp
*
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/flow/ShortestSuccessivePath.hpp>
namespace NetworKit {
MinFlowShortestSuccessivePath::MinFlowShortestSuccessivePath(
    const Graph &G,
    const std::string &capacityName,
    const std::string &supplyName)
    : graph(G), capacityAttributeName(capacityName), supplyAttributeName(supplyName) {

    if (!G.isDirected()) {
        throw std::runtime_error(
            "MinFlowShortestSuccessivePath: Graph must be directed");
    }

    if (!G.isWeighted()) {
        throw std::runtime_error(
            "MinFlowShortestSuccessivePath: Graph must be weighted.");
    }

    try {
        // this internally calls AttributeMap::find and will throw if the attribute does not exist
       (void)graph.edgeAttributes().find(capacityName);
    } catch (const std::runtime_error &e) {
        // catch the generic “No such attribute” and give our own
        throw std::runtime_error(
            "MinFlowShortestSuccessivePath: Provided edge attribute '" + capacityName + "' not found");
    }

    try {
        // this internally calls AttributeMap::find and will throw if the attribute does not exist
        (void)graph.nodeAttributes().find(supplyName);
    } catch (const std::runtime_error &e) {
        // catch the generic “No such attribute” and give our own
        throw std::runtime_error(
            "MinFlowShortestSuccessivePath: Provided node attribute '" + supplyName + "' not found");
    }
}


}