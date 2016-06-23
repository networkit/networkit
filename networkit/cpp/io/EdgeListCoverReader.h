#ifndef EDGELISTCOVERREADER_H
#define EDGELISTCOVERREADER_H

#include "CoverReader.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class EdgeListCoverReader : public CoverReader {
	public:
		/**
		 * Constructs the EdgeListCoverReader class with @a firstNode as the index of the first node in the file.
		 * @param[in]	firstNode	Index of the first node in the file.
		 */
		EdgeListCoverReader(node firstNode = 1);

		/**
		 * Read a cover from a file. File format: each line contains the node ids of one subset.
		 *
		 * @param[in]	path	The path to the input file
		 * @param[in]	G		The graph for which the cover shall be read
		 * @return The cover instance
		 */
		virtual Cover read(std::string path, Graph& G);
	private:
		node firstNode;
};
}

#endif // EDGELISTCOVERREADER_H
