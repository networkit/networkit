#ifndef NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_

#include <concepts>
#include <span>

#include <networkit/Globals.hpp>

namespace NetworKit::IsomorphismDetails {

/// Receives a complete mapping, indexed by pattern node, from a search implementation. Returns
/// false once the search must stop because the cap on the number of matches is reached. The search
/// implementations take the reporter as a template parameter, so that the call per match needs no
/// type erasure.
template <typename Reporter>
concept MatchReporter = std::predicate<Reporter &, std::span<const node>>;

} // namespace NetworKit::IsomorphismDetails

#endif // NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
