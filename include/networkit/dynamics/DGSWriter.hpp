/*
 * DGSWriter.h
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef NETWORKIT_DYNAMICS_DGS_WRITER_HPP_
#define NETWORKIT_DYNAMICS_DGS_WRITER_HPP_

#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class DGSWriter {
public:
    DGSWriter() = default;

    void write(std::vector<GraphEvent>& stream, std::string path);

};

} /* namespace NetworKit */

#endif // NETWORKIT_DYNAMICS_DGS_WRITER_HPP_
