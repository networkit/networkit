/*
 * ImageGraphWriter.h
 *
 *  Created on: 26.04.2014
 *      Author: Michael
 */

#ifndef IMAGEGRAPHWRITER_H_
#define IMAGEGRAPHWRITER_H_

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include <vector>

#include "../imageio/png++/image.hpp"
#include <sstream>

namespace NetworKit {

/**
 * The ImageGraphWriter class writes a partitioned Graph (constructed from a png image) to disk
 */
class ImageGraphWriter {
public:
	/** Default constructor */
	ImageGraphWriter();
	/** Destructor */
	virtual ~ImageGraphWriter();

	/**
	 * Writes the partitions in @a partition into one single image - each partition gets painted with the mean color within that partition.
	 * The @a partition belongs to @a graph which has been constructed from a png image located at @a pathToSourceImage.
	 * @param graph The graph.
	 * @param partition The partition/clustering of @a graph.
	 * @param pathToSourceImage The path to the image which the @a graph represents.
	 */
	void writePartitionsToSingleImage(Graph &graph, const Partition &partition, const std::string &pathToSourceImage) const;

	/**
	 * Writes each partition in @a partition into an image file. The partition part is copied from the source image.
	 * @param graph The graph.
	 * @param partition The partition/clustering of @a graph.
	 * @param pathToSourceImage The path to the image which the @a graph represents.
	 */
	void writePartitions(Graph &graph, const Partition &partition, const std::string &pathToSourceImage) const;

private:
	/**
	 * Creates an image from a partition of @a graph specified by @a partitionMembers.
	 * @param graph The graph.
	 * @param partitionMembers The elements of the partition.
	 * @param sourceImage The source image.
	 */
	png::image<png::rgb_pixel> createImageForPartition(Graph &graph, const std::set<index> partitionMembers, const png::image<png::rgb_pixel> &sourceImage) const;
};

} /* namespace NetworKit */

#endif /* IMAGEGRAPHWRITER_H_ */
