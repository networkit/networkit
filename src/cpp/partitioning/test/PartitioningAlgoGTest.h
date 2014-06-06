#ifndef PARTITIONINGALGOGTEST_H_
#define PARTITIONINGALGOGTEST_H_

#ifndef NOGTEST

#include <gtest/gtest.h>

#include "../BalancedLabelPropagation.h"
#include "../../auxiliary/Log.h"
#include "../../community/Modularity.h"
#include "../../community/EdgeCut.h"
#include "../../graph/GraphGenerator.h"
#include "../../community/ClusteringGenerator.h"
#include "../../structures/Partition.h"
#include "../../io/METISGraphReader.h"
#include "../../io/DibapGraphReader.h"
#include "../../viz/PostscriptWriter.h"

namespace NetworKit {

class PartitioningAlgoGTest: public testing::Test {
};

} /* namespace NetworKit */

#endif /* PARTITIONINGALGOGTEST_H_ */

#endif /*NOGTEST */
