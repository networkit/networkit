/*
 * Clustering.h
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

namespace EnsembleClustering {

// TODO: import
typedef int Cluster;
typedef int Node;

class Clustering {

public:

	Clustering();

	virtual ~Clustering();

	virtual Cluster getCluster(Node u);

	virtual void addToCluster(Cluster c, Node u);

	virtual void addToNewCluster(Node u);

	virtual void moveToCluster(Cluster c, Node u);



};

} /* namespace EnsembleClustering */
#endif /* CLUSTERING_H_ */
