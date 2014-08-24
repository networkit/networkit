#ifndef COVERREADER_H_
#define COVERREADER_H_

#include "../structures/Cover.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class CoverReader {
	public:
		/**
		 * Read a cover from a file. File format: each line contains the node ids of one subset.
		 *
		 * @param[in]	path	The path to the input file
		 * @param[in]	G		The graph for which the cover shall be read
		 * @return The cover instance
		 */
		virtual Cover read(std::string path, Graph& G);

};

} /* namespace NetworKit */

#endif // COVERREADER_H_
