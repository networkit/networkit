#ifndef NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_

// Private header of the isomorphism module. Not installed, not part of the public API.

#include <functional>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

namespace NetworKit {
namespace IsomorphismDetails {

/**
 * Callback through which a search implementation (VF2Impl, VF3Impl, RIImpl) reports a complete
 * mapping, indexed by pattern node. It returns false once the search must stop because the cap on
 * the number of matches is reached.
 *
 * The implementations are not subclasses of SubgraphIsomorphism, so run() passes a lambda that
 * calls the protected SubgraphIsomorphism::reportMatch():
 *
 * @code
 * VF2Impl(..., [this](const SubgraphIsomorphism::Match &m) { return reportMatch(m); }).run();
 * @endcode
 */
using MatchReporter = std::function<bool(const SubgraphIsomorphism::Match &)>;

} // namespace IsomorphismDetails
} // namespace NetworKit

#endif // NETWORKIT_CPP_ISOMORPHISM_MATCH_REPORTER_HPP_
