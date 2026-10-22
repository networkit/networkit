#include <stdexcept>
#include <utility>

#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {

SubgraphIsomorphism::SubgraphIsomorphism(const Graph &pattern, const Graph &target,
                                         Semantics semantics, count maxMatches)
    : Algorithm(), pattern(&pattern), target(&target), semantics(semantics), maxMatches(maxMatches),
      matchCount(0), storeMatches(true), matchesStored(false) {
    // Python passes the semantics as an integer, so any value can arrive here.
    if (semantics != Semantics::INDUCED && semantics != Semantics::MONOMORPHISM)
        throw std::invalid_argument("Unknown subgraph isomorphism semantics");

    validateInput();
}

void SubgraphIsomorphism::validateInput() const {
    if (pattern->isDirected() != target->isDirected())
        throw std::runtime_error(
            "Pattern and target graph must either both be directed or both be undirected");

    // Target self-loops are harmless, since both semantics only constrain pairs of distinct nodes.
    if (pattern->numberOfSelfLoops() > 0)
        throw std::runtime_error("Subgraph isomorphism is undefined for patterns with self-loops");
}

void SubgraphIsomorphism::validateNodeLabels(const std::vector<index> &patternNodeLabels,
                                             const std::vector<index> &targetNodeLabels) const {
    // Two empty vectors clear the labels.
    if (patternNodeLabels.empty() && targetNodeLabels.empty())
        return;

    if (patternNodeLabels.size() < pattern->upperNodeIdBound())
        throw std::runtime_error("Pattern label vector is shorter than the pattern's "
                                 "upperNodeIdBound()");

    if (targetNodeLabels.size() < target->upperNodeIdBound())
        throw std::runtime_error("Target label vector is shorter than the target's "
                                 "upperNodeIdBound()");
}

void SubgraphIsomorphism::validateEdgeLabels(const std::vector<index> &patternEdgeLabels,
                                             const std::vector<index> &targetEdgeLabels) const {
    if (patternEdgeLabels.empty() && targetEdgeLabels.empty())
        return;

    if (!pattern->hasEdgeIds())
        throw std::runtime_error("Pattern graph has no edge ids - call indexEdges() on it before "
                                 "setting edge labels");

    if (!target->hasEdgeIds())
        throw std::runtime_error("Target graph has no edge ids - call indexEdges() on it before "
                                 "setting edge labels");

    if (patternEdgeLabels.size() < pattern->upperEdgeIdBound())
        throw std::runtime_error("Pattern edge label vector is shorter than the pattern's "
                                 "upperEdgeIdBound()");

    if (targetEdgeLabels.size() < target->upperEdgeIdBound())
        throw std::runtime_error("Target edge label vector is shorter than the target's "
                                 "upperEdgeIdBound()");
}

void SubgraphIsomorphism::setNodeLabels(const std::vector<index> &patternNodeLabels,
                                        const std::vector<index> &targetNodeLabels) {
    // Validate first, so that a rejected call changes nothing.
    validateNodeLabels(patternNodeLabels, targetNodeLabels);

    this->patternNodeLabels = patternNodeLabels;
    this->targetNodeLabels = targetNodeLabels;
}

void SubgraphIsomorphism::setEdgeLabels(const std::vector<index> &patternEdgeLabels,
                                        const std::vector<index> &targetEdgeLabels) {
    validateEdgeLabels(patternEdgeLabels, targetEdgeLabels);

    this->patternEdgeLabels = patternEdgeLabels;
    this->targetEdgeLabels = targetEdgeLabels;
}

void SubgraphIsomorphism::setCallback(MatchCallback callback) {
    this->callback = std::move(callback);
}

void SubgraphIsomorphism::setStoreMatches(bool storeMatches) {
    this->storeMatches = storeMatches;
}

void SubgraphIsomorphism::prepareRun() {
    hasRun = false;

    result.clear();
    result.shrink_to_fit();

    matchCount = 0;
    matchesStored = storeMatches && !hasCallback();

    // The graphs may have changed since construction, which can leave the label vectors too
    // short. The check runs after the reset, so that a rejected run leaves no results behind.
    validateInput();
    validateNodeLabels(patternNodeLabels, targetNodeLabels);
    validateEdgeLabels(patternEdgeLabels, targetEdgeLabels);
}

bool SubgraphIsomorphism::reportMatch(const Match &match) {
    ++matchCount;

    if (callback)
        callback(match);
    else if (matchesStored)
        result.push_back(match);

    return maxMatches == 0 || matchCount < maxMatches;
}

void SubgraphIsomorphism::finishRun() {
    // reportMatch() already enforced the cap.
    hasRun = true;
}

const std::vector<SubgraphIsomorphism::Match> &SubgraphIsomorphism::getMatches() const {
    assureFinished();
    if (!matchesStored)
        throw std::runtime_error("SubgraphIsomorphism did not store the matches of the last run, "
                                 "since a callback was set or setStoreMatches(false) was called");
    return result;
}

count SubgraphIsomorphism::numberOfMatches() const {
    assureFinished();
    return matchCount;
}

bool SubgraphIsomorphism::hasMatch() const {
    assureFinished();
    return matchCount > 0;
}

} // namespace NetworKit
