#ifndef CLUSTERPRODUCT_H
#define CLUSTERPRODUCT_H

#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Class for calculating the intersection of two partitions, i.e. the clustering with the fewest clusters
 * such that each cluster is a subset of a cluster in both partitions.
 */
class PartitionIntersection {
	public:
		/**
		 * Calculate the intersection of two partitions @a zeta and @a eta
		 * @param zeta	The first partition
		 * @param eta	The second partition
		 * @return The intersection of @a zeta and @a eta
		 */
		Partition calculate(const Partition &zeta, const Partition &eta);
};

} /* namespace NetworKit */


#endif // CLUSTERPRODUCT_H
