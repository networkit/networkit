/*
 * DynDynWeaklyConnectedComponents.cpp
 *
 *  Created on: June 20, 2017
 *      Author: Eugenio Angriman
 */

#include "DynConnectedComponentsImpl.hpp"

#include <span>

#include <networkit/components/DynWeaklyConnectedComponents.hpp>

namespace NetworKit {
DynWeaklyConnectedComponents::DynWeaklyConnectedComponents(const Graph &G)
    : ComponentDecomposition(G),
      impl(new DynConnectedComponentsDetails::DynConnectedComponentsImpl<true>{G, component}) {}

DynWeaklyConnectedComponents::~DynWeaklyConnectedComponents() = default;

void DynWeaklyConnectedComponents::run() {
    impl->run();
    hasRun = true;
}

void DynWeaklyConnectedComponents::update(GraphEvent event) {
    impl->update(event);
}

void DynWeaklyConnectedComponents::updateBatch(std::span<const GraphEvent> batch) {
    impl->updateBatch(batch);
}
} // namespace NetworKit
