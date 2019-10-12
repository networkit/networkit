/*
 * PostscriptWriter.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */
// networkit-format

#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/viz/PostscriptWriter.hpp>

namespace NetworKit {

namespace PostscriptWriterColors {
struct RGBColor {
    float r;
    float g;
    float b;
};

static RGBColor fromCyclicRotation(size_t index) {
    constexpr RGBColor colors[] = {
        RGBColor({1.0f, 0.0f, 0.0f}), RGBColor({1.0f, 0.5f, 0.0f}), RGBColor({1.0f, 1.0f, 0.0f}),
        RGBColor({0.5f, 1.0f, 0.0f}), RGBColor({0.0f, 1.0f, 0.0f}), RGBColor({0.0f, 1.0f, 0.5f}),
        RGBColor({0.0f, 1.0f, 1.0f}), RGBColor({0.0f, 0.5f, 1.0f}), RGBColor({0.0f, 0.0f, 1.0f}),
        RGBColor({0.5f, 0.0f, 1.0f}), RGBColor({1.0f, 0.0f, 1.0f}), RGBColor({1.0f, 0.0f, 0.5f}),
        RGBColor({0.6f, 0.0f, 0.0f}), RGBColor({0.6f, 0.3f, 0.0f}), RGBColor({0.6f, 0.6f, 0.0f}),
        RGBColor({0.3f, 0.6f, 0.0f}), RGBColor({0.0f, 0.6f, 0.0f}), RGBColor({0.0f, 0.6f, 0.3f}),
        RGBColor({0.0f, 0.6f, 0.6f}), RGBColor({0.0f, 0.3f, 0.6f}), RGBColor({0.0f, 0.0f, 0.6f}),
        RGBColor({0.3f, 0.0f, 0.6f}), RGBColor({0.6f, 0.0f, 0.6f}), RGBColor({0.6f, 0.0f, 0.3f})};
    constexpr auto numColors = sizeof(colors) / sizeof(colors[0]);

    return colors[index % numColors];
}
} // namespace PostscriptWriterColors

PostscriptWriter::PostscriptWriter(bool isTorus) : wrapAround(isTorus), ps_size{1020.0, 1020.0} {}

void PostscriptWriter::computeBoundaryBox(const std::vector<coord2d> &coordinates) {
    ps_min = {std::numeric_limits<coord>::max(), std::numeric_limits<coord>::max()};
    ps_max = {std::numeric_limits<coord>::min(), std::numeric_limits<coord>::min()};
    for (const auto p : coordinates) {
        ps_min.first = std::min(ps_min.first, p.first);
        ps_max.first = std::max(ps_max.first, p.first);
        ps_min.second = std::min(ps_min.second, p.second);
        ps_max.second = std::max(ps_max.second, p.second);
    }

    ps_scale = {(ps_size.first - 2 * ps_border.first) / (ps_max.first - ps_min.first),
                (ps_size.second - 2 * ps_border.second) / (ps_max.second - ps_min.second)};
}

void PostscriptWriter::writeHeader(std::ofstream &file) const {
    /* Header */
    if (wrapAround) {
        file << "%!PS-Adobe-3.0 EPSF-3.0\n";
    } else {
        file << "%!PS-Adobe-1.0\n";
    }
    file << "%%Title: NetworKit visualization\n";
    file << "%%BoundingBox: 0.000 0.000 " << ps_size.first << " " << ps_size.second << "\n";
    file << "%%EndComments\n";
    if (!wrapAround) {
        file << "%%EndProlog\n";
        file << "gsave\n";
    }
}

void PostscriptWriter::writeMacros(std::ofstream &file) const {
    /* Macros */
    file << "/p {newpath} bind def\n"
            "/m {moveto} bind def\n"
            "/r {rmoveto} bind def\n"
            "/k {rlineto} bind def\n"
            "/l {lineto} bind def\n"
            "/n {rlineto} bind def\n"
            "/c {setrgbcolor} bind def\n"
            "/s {stroke} bind def\n"
            "/w {setlinewidth} bind def\n"
            "/h {show} bind def\n"
            "/a {arc closepath fill} bind def\n"
            "/b {closepath eofill} bind def\n";
}

void PostscriptWriter::writeClustering(const Graph &g, const std::vector<coord2d> &coordinates,
                                       const Partition &clustering, std::ofstream &file) {

    auto adjustToBoundingBox([&](coord2d p) {
        return coord2d{(p.first - ps_min.first) * ps_scale.first + ps_border.first,
                       (p.second - ps_min.second) * ps_scale.second + ps_border.second};
    });

    auto adjustWrapAround = [](coord2d p) {
        auto adjust = [](coord val) { // TODO: externalize constants
            if (val > 500.0f)
                return val - 1000.0f;
            if (val < -500.0f)
                return val + 1000.0f;
            return val;
        };

        return coord2d{adjust(p.first), adjust(p.second)};
    };

    g.forEdges([&](node u, node v) {
        // set edge color
        if (clustering[u] == clustering[v] && clustering[u] != none) {
            // same cluster
            const auto color = PostscriptWriterColors::fromCyclicRotation(clustering[u]);
            file << color.r << " " << color.g << " " << color.b << " c ";
        } else {
            // different clusters -> grey
            file << "0.80 0.80 0.80 c 1.0 w ";
        }

        // set edge start and end point
        auto start = adjustToBoundingBox(coordinates[u]);
        auto end = adjustToBoundingBox(coordinates[v]);
        auto diff = coord2d{end.first - start.first, end.second - start.second};

        if (wrapAround) {
            diff = adjustWrapAround(diff);
            end = {start.first + diff.first, start.second + diff.second};
        }

        // write edge to file
        file << "p " << start.first << " " << start.second << " m " << end.first << " "
             << end.second << " l s\n";
    });

    // draw vertices
    float dotsize = 2.0;
    g.forNodes([&](node u) {
        if (clustering[u] != none) {
            // change color
            const auto color = PostscriptWriterColors::fromCyclicRotation(clustering[u]);
            file << color.r << " " << color.g << " " << color.b << " c ";
        } else {
            file << "0.0 0.0 0.0 c ";
        }

        const auto point = adjustToBoundingBox(coordinates[u]);

        file << "p " << point.first << " " << point.second << " " << dotsize
             << " 0.00 360.00 a s\n";
    });
}

void PostscriptWriter::init(std::ofstream &file) const {
    file.precision(3);
    file << std::fixed;

    writeHeader(file);
    writeMacros(file);
    file << "0.000 0.000 0.000 c\n";
}

void PostscriptWriter::write(const Graph &g, const std::vector<coord2d> &coordinates,
                             const Partition &clustering, const std::string &path) {
    std::ofstream file(path);
    init(file);
    computeBoundaryBox(coordinates);

    writeClustering(g, coordinates, clustering, file);

    if (!wrapAround) {
        file << "grestore\n";
    }

    file.close();
}

void PostscriptWriter::write(const Graph &g, const std::vector<coord2d> &coordinates,
                             const std::string &path) {
    ClusteringGenerator gen;
    Partition allNone = gen.makeOneClustering(g);
    write(g, coordinates, allNone, path);
}

} /* namespace NetworKit */
