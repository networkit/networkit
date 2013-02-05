/*
 * ClusteringWriter.h
 *
 *  Created on: 22.01.2013
 *      Author: cls
 */

#ifndef CLUSTERINGWRITER_H_
#define CLUSTERINGWRITER_H_

#include <fstream>

#include "../clustering/base/Clustering.h"

namespace EnsembleClustering {

class ClusteringWriter {

public:

	ClusteringWriter();

	virtual ~ClusteringWriter();

	virtual void write(Clustering& zeta, std::string path);
};

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGWRITER_H_ */
