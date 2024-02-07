#ifndef NETWORKIT_GENERATORS_CONFIGURATION_MODEL_HPP_
#define NETWORKIT_GENERATORS_CONFIGURATION_MODEL_HPP_

#include <set>
#include <vector>
#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Graph generator for generating a random simple graph with exactly the given degree sequence.
 */
class ConfigurationModel final : public StaticDegreeSequenceGenerator {
public:
    /**
     * Initializes the configuration model generator with the given degree sequence @a sequence.
     *
     * @param sequence std::vector<count> with the input sequence of degrees.
     */
    ConfigurationModel(const std::vector<count> &sequence);

    /**
     * Generates a graph according to the configuration model.
     *
     * @throws std::runtime_error If the sequence cannot be generated.
     * @return The generated graph
     */
    Graph generate() override;
};

} // namespace NetworKit

#endif // NETWORKIT_GENERATORS_CONFIGURATION_MODEL_HPP_
