#ifndef NETWORKIT_IO_COVER_WRITER_HPP_
#define NETWORKIT_IO_COVER_WRITER_HPP_

#include <networkit/structures/Cover.hpp>
#include <networkit/graph/Graph.hpp>

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

#endif // NETWORKIT_IO_COVER_WRITER_HPP_
