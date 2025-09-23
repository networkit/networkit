#ifndef NETWORKIT_BASE_DYN_ALGORITHM_HPP_
#define NETWORKIT_BASE_DYN_ALGORITHM_HPP_

#include <span>

#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

class DynAlgorithm {

public:
    /**
     * Virtual default destructor
     */
    virtual ~DynAlgorithm() = default;

    /**
     * The generic update method for updating data structure after an update.
     */
    virtual void update(GraphEvent e) = 0;

    /**
     * The generic update method for updating data structure after a batch of updates.
     */
    virtual void updateBatch(std::span<const GraphEvent> batch) = 0;
};

} // namespace NetworKit

#endif // NETWORKIT_BASE_DYN_ALGORITHM_HPP_
