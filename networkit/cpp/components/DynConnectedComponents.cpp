/*
 * DynConnectedComponents.cpp
 *
 *  Created on: June 2017
 *      Author: Eugenio Angriman
 */

#include "DynConnectedComponentsImpl.hpp"

#include <networkit/components/DynConnectedComponents.hpp>

namespace NetworKit {

DynConnectedComponents::DynConnectedComponents(const Graph &G)
    : ComponentDecomposition(G),
      impl(new DynConnectedComponentsDetails::DynConnectedComponentsImpl<false>{G, component}) {}

DynConnectedComponents::~DynConnectedComponents() = default;

void DynConnectedComponents::run() {
    impl->run();
    hasRun = true;
}

void DynConnectedComponents::update(GraphEvent event) {
    impl->update(event);
}

void DynConnectedComponents::updateBatch(const std::vector<GraphEvent> &batch) {
    impl->updateBatch(batch);
}
} // namespace NetworKit
