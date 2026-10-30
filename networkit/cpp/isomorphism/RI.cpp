#include <span>
#include <stdexcept>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/isomorphism/RI.hpp>

#include "RIImpl.hpp"

namespace NetworKit {

RI::RI(const Graph &pattern, const Graph &target, Variant variant, Semantics semantics,
       count maxMatches, double domainFilterThreshold)
    : SubgraphIsomorphism(pattern, target, semantics, maxMatches), variant(variant),
      domainFilterThreshold(domainFilterThreshold) {
    // Python passes the variant as an integer, so any value can arrive here.
    if (variant != Variant::RI && variant != Variant::RI_DS)
        throw std::invalid_argument("Unknown RI variant");

    // Negated, so that NaN fails too.
    if (!(domainFilterThreshold >= 0.0 && domainFilterThreshold <= 1.0))
        throw std::invalid_argument("The domain filter threshold must be between 0 and 1");
}

void RI::run() {
    using IsomorphismDetails::prepareRISearch;
    using IsomorphismDetails::RIImpl;
    using IsomorphismDetails::RISearchSetup;

    Aux::SignalHandler handler;
    prepareRun();

    const RISearchSetup setup =
        prepareRISearch(*pattern, *target, patternNodeLabels, targetNodeLabels, patternEdgeLabels,
                        targetEdgeLabels, variant, domainFilterThreshold, "RI");

    RIImpl(setup.patternGraph, setup.targetGraph, patternNodeLabels, targetNodeLabels,
           setup.ordering, setup.domains, semantics, handler,
           [this](std::span<const node> match) { return reportMatch(match); })
        .run();

    // RIImpl only polls isRunning(), so an interrupted search returns normally and throws here.
    handler.assureRunning();
    finishRun();
}

} // namespace NetworKit
