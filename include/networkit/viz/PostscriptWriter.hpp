/*
 * PostscriptWriter.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */
// networkit-format

#ifndef NETWORKIT_VIZ_POSTSCRIPT_WRITER_HPP_
#define NETWORKIT_VIZ_POSTSCRIPT_WRITER_HPP_

#include <fstream>
#include <string>
#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup viz
 * EPS output of graphs with 2D coordinates
 */
class PostscriptWriter {
public:
    /**
     * @param[in] isTorus Specifies whether the visualization square is treated as torus,
     * i.e. with wrap-around boundaries (edge can leave the square and enter at the opposite
     * side. By default, it is set to false.
     */
    PostscriptWriter(bool isTorus = false);

    /**
     * Outputs an EPS file with name @a filename of the graph @a g with 2D coordinates.
     * The colors are chosen to visualize the specified @a clustering.
     * @param[in] g Graph to be visualized.
     * @param[in] clustering Clustering of the graph, visualized by different colors.
     * @param[in] filename Name of file to write to.
     */
    void write(const Graph &g, const std::vector<coord2d> &coordinates, const Partition &clustering,
               const std::string &filename);

    /**
     * Outputs an EPS file with name @a filename of the graph @a g with 2D coordinates.
     * @param[in] g Graph to be visualized.
     * @param[in] filename Name of file to write to.
     */
    void write(const Graph &g, const std::vector<coord2d> &coordinates,
               const std::string &filename);

private:
    bool wrapAround;

    coord2d ps_size;
    coord2d ps_border;

    coord2d ps_min;
    coord2d ps_max;
    coord2d ps_scale;

    void init(std::ofstream &file) const;

    void computeBoundaryBox(const std::vector<coord2d> &coordinates);

    void writeHeader(std::ofstream &file) const;
    void writeMacros(std::ofstream &file) const;
    void writeClustering(const Graph &g, const std::vector<coord2d> &coordinates,
                         const Partition &clustering, std::ofstream &file);
};

} /* namespace NetworKit */
#endif // NETWORKIT_VIZ_POSTSCRIPT_WRITER_HPP_
