#ifndef NETWORKIT_COMMUNITY_PARTITION_FRAGMENTATION_HPP_
#define NETWORKIT_COMMUNITY_PARTITION_FRAGMENTATION_HPP_

#include <networkit/community/LocalPartitionEvaluation.hpp>

namespace NetworKit {

/**
 * This measure evaluates how fragmented a partition is. The fragmentation of a single cluster is defined as one minus the
 * number of nodes in its maximum connected componented divided by its total number of nodes. Smaller values thus indicate a smaller fragmentation.
 */
class PartitionFragmentation : public LocalPartitionEvaluation {
public:
    using LocalPartitionEvaluation::LocalPartitionEvaluation;

    /**
     * Execute the algorithm.
     */
    virtual void run() override;

    /**
     * @return true - smaller values are better than larger values.
     */
    virtual bool isSmallBetter() const override { return true; };

    /**
     * @return false - only minor parts of this implementation are parallel.
     */
    virtual bool isParallel() const override { return false; };
};

}

#endif // NETWORKIT_COMMUNITY_PARTITION_FRAGMENTATION_HPP_
