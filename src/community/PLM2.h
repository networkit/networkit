/*
 * PLM2.h
 *
 *  Created on: 22.10.2013
 *      Author: cls
 */

 


#ifndef PLM2_H_
#define PLM2_H_

#include "Clusterer.h"



namespace NetworKit {

/**
 * Lock-free parallel implementation of the Louvain method
 * for community detection in networks.
 */
class PLM2: public NetworKit::Clusterer {
public:

	/**
	 * @param[in]	par		parallelization strategy
	 * @param[in]	gamma	multi-resolution modularity parameter:
	 * 							1.0 -> standard modularity
	 * 							0.0 -> one community
	 * 							2m 	-> singleton communities
	 *
	 */
	PLM2(std::string par="simple", double gamma = 1.0);

	virtual ~PLM2();

	/**
	 * Perform one optimization pass.
	 */
	virtual Clustering pass(const Graph& G);


	/**
	 * Run the algorithm.
	 */
	Clustering run(Graph& G) override;


	/**
	 * @return string representation of algorithm and parameters.
	 */
	std::string toString() const override;

protected:

	std::string parallelism; //!< switch for the kind of parallelization strategy to use
	double gamma;	//!< multi-resolution modularity parameter
	bool anyChange;	//!< indicates whether any change was made to the clustering in the last pass over the nodes

};

} /* namespace NetworKit */
#endif /* PLM2_H_ */


