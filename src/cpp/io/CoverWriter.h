#ifndef COVERWRITER_H
#define COVERWRITER_H

#include "../structures/Cover.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup io
 * Write a clustering to a file.
 */

class CoverWriter
{
	public:

		virtual void write(Cover& zeta, const std::string& path) const;
};
}

#endif // COVERWRITER_H
