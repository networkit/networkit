#ifndef NETWORKIT_GENERATORS_CONFIGURATION_MODEL_HPP_
#define NETWORKIT_GENERATORS_CONFIGURATION_MODEL_HPP_

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
     * @param ignoreIfNotRealizable Ignore if the sequence is not realizable and generate a graph
     *        anyway.
     */
    ConfigurationModel(const std::vector<count> &sequence);

    /**
     * Generate a graph according to the configuration model.
     *
     * Issues a INFO log message if the wanted number of edge swaps cannot be performed because of
     * the limit of attempts (see in the description of the class for details).
     *
     * @throws std::runtime_error If the sequence cannot be generated.
     * @return The generated graph
     */
    Graph generate() override;
};

} // namespace NetworKit

#endif // NETWORKIT_GENERATORS_CONFIGURATION_MODEL_HPP_
