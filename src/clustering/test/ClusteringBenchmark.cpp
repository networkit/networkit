/*
 * ClusteringBenchmark.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "ClusteringBenchmark.h"

#include "../Clustering.h"
#include "../../structures/Partition.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Random.h"


namespace NetworKit {

TEST_F(ClusteringBenchmark, PartitionVsClustering) {
	count n = 1e5;
	n *= 2;
	count merge_ops = n/2;
	std::vector<index> move_ids;
	std::vector<index> merge_ids;
	for (index i = 0; i < merge_ops; i++) {
		merge_ids.push_back(Aux::Random::integer(n-1));
		move_ids.push_back(Aux::Random::integer(n-1));		
	}
	
	Aux::Timer pTimer;
	pTimer.start();
	Partition p(n);
	p.allToSingletons();
	for (index i = 0; i < merge_ops; ++i) {
		p.mergeSubsets(p[merge_ids[i]],p[move_ids[i]]);
	}
	p.subsetSizes();
	p.compact();
	for (index i = 0; i < merge_ops; ++i) {
		p.moveToSubset(p[merge_ids[i]],move_ids[i]);
	}
	p.subsetSizeMap();
	pTimer.stop();
	
	Aux::Timer cTimer;
	cTimer.start();
	Clustering c(n);
	c.allToSingletons();
	for (index i = 0; i < merge_ops; ++i) {
		c.mergeClusters(c[merge_ids[i]],c[move_ids[i]]);
	}
	c.clusterSizes();
	c.compact();
	for (index i = 0; i < merge_ops; ++i) {
		c.moveToCluster(c[merge_ids[i]],move_ids[i]);
	}
	c.clusterSizeMap();
	cTimer.stop();

	std::cout<<"Benchmarking Partition and Clustering with:"<<std::endl;
	std::cout<<"\t"<<n<<" elements"<<std::endl;
	std::cout<<"\t"<<merge_ops<<" merge operations"<<std::endl;
	std::cout<<"\tclusterSizes()"<<std::endl;
	std::cout<<"\tcompaction"<<std::endl;
	std::cout<<"\t"<<merge_ops<<" move operations"<<std::endl;
	std::cout<<"\tclusterSizeMap()"<<std::endl;
	std::cout<<"Time for Partition:\t"<<pTimer.elapsedTag()<<std::endl;
	std::cout<<"Time for Clustering:\t"<<cTimer.elapsedTag()<<std::endl;
	
}} /* namespace NetworKit */

#endif /*NOGTEST */

