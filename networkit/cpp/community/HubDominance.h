#ifndef HUBDOMINANCE_H
#define HUBDOMINANCE_H

#include "../structures/Partition.h"
#include "../structures/Cover.h"
#include "QualityMeasure.h"

namespace NetworKit {

/**
 * A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
 * cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
 * the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
 * value for all clusters is defined as the average of all clusters.
 * Strictly speaking this is not a quality measure as this is rather dependent on the type of the
 * considered graph, for more information see
 * Lancichinetti A, Kivelä M, Saramäki J, Fortunato S (2010)
 * Characterizing the Community Structure of Complex Networks
 * PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
 * http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976
 */
class HubDominance : public QualityMeasure {
public:
	/**
	 * Calculates the dominance of hubs in the given Partition @a zeta of the given
	 * Graph @a G.
	 *
	 * @param zeta	The partition for which the hub dominance shall be calculated
	 * @param G	The graph that is partitioned in @a zeta
	 * @return The average hub dominance of @a zeta
	 */
	virtual double getQuality(const Partition& zeta, const Graph& G) override;
	/**
	 * Calculates the dominance of hubs in the given Cover @a zeta of the given
	 * Graph @a G.
	 *
	 * @param zeta	The cover for which the hub dominance shall be calculated
	 * @param G	The graph that is partitioned in @a zeta
	 * @return The average hub dominance of @a zeta
	 */
	virtual double getQuality(const Cover& zeta, const Graph& G);
};

} // namespace NetworKit

#endif // HUBDOMINANCE_H
