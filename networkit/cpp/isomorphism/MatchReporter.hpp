#ifndef NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_

#include <functional>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {
namespace IsomorphismDetails {

/// Receives a complete mapping, indexed by pattern node, from a search implementation. Returns
/// false once the search must stop because the cap on the number of matches is reached.
using MatchReporter = std::function<bool(const SubgraphIsomorphism::Match &)>;

} // namespace IsomorphismDetails
} // namespace NetworKit

#endif // NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
