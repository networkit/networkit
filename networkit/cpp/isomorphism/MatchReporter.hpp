#ifndef NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_

#include <functional>
#include <span>

#include <networkit/Globals.hpp>

namespace NetworKit {
namespace IsomorphismDetails {

/// Receives a complete mapping, indexed by pattern node, from a search implementation. Returns
/// false once the search must stop because the cap on the number of matches is reached.
using MatchReporter = std::function<bool(std::span<const node>)>;

} // namespace IsomorphismDetails
} // namespace NetworKit

#endif // NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
