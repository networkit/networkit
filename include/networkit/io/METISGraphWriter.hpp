/*
 * METISGraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef METISGRAPHWRITER_H_
#define METISGRAPHWRITER_H_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class METISGraphWriter final : public GraphWriter {

public:

	METISGraphWriter() = default;

	void write(const Graph &G, const std::string &path) override;

	void write(const Graph &G, bool weighted, const std::string &path);
};

} /* namespace NetworKit */
#endif /* METISGRAPHWRITER_H_ */
