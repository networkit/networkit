/*
 * KernighanLin.h
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#ifndef KERNIGHANLIN_H_
#define KERNIGHANLIN_H_

#include "BalancedPartitioner.h"

namespace NetworKit {

class KernighanLin: public BalancedPartitioner {
public:
	KernighanLin();
	virtual ~KernighanLin();

	virtual Clustering run(Graph& G, count numBlocks);
	virtual Clustering& rerun(Graph& G, count numBlocks, Clustering& partition);

	virtual Clustering& postsmooth(Graph& graph, count numBlocks, Clustering& partition);
};

} /* namespace NetworKit */
#endif /* KERNIGHANLIN_H_ */
