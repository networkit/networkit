/*
 * ImageSegmentation.h
 *
 *  Created on: 27.04.2014
 *      Author: Michael
 */

#ifndef IMAGESEGMENTATION_H_
#define IMAGESEGMENTATION_H_

#include "../community/CommunityDetectionAlgorithm.h"
#include "ImageGraphReader.h"
#include "ImageGraphWriter.h"

namespace NetworKit {

/**
 * The ImageSegmentation class bundles the steps of segmenting an image.
 */
class ImageSegmentation {
public:

	/** Default constructor */
	ImageSegmentation();
	/** Destructor */
	virtual ~ImageSegmentation();

	/**
	 * Creates a segmentation of the image located at @a path. The clustering algorithm is specified by @a algorithm.
	 * @param algorithm The algorithm used for clustering.
	 * @param path The path to the image which should be segmented.
	 * @param allParititonsToOneFile If true, all found partitions will be written to one file rather than each partition in an extra image file.
	 */
	void createSegmentation(CommunityDetectionAlgorithm &algorithm, const std::string &path, bool allPartitionsToOneFile = true);
};

} /* namespace NetworKit */

#endif /* IMAGESEGMENTATION_H_ */
