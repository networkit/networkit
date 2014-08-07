#ifndef CLUSTERPRODUCT_H
#define CLUSTERPRODUCT_H

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * The product of two partitions is defined as the partitions where each cluster is the intersection
 * of a cluster in the first and in the second clustering
 */
// TODO: rename to PartitionProduct
class ClusteringProduct {
	public:
		/**
		 * Calculate the product of two partitions @a zeta and @a eta
		 * @param zeta	The first partition
		 * @param eta	The second partition
		 * @return The product of @a zeta and @a eta
		 */
		Partition calculate(const Partition &zeta, const Partition &eta);
};

} /* namespace NetworKit */


#endif // CLUSTERPRODUCT_H
