/*
 * NetworkitBinaryWriterImpl.hpp
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_WRITER_IMPL_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_WRITER_IMPL_HPP_

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <unordered_map>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>

namespace NetworKit {

template <class GraphT>
void NetworkitBinaryWriter::write(const GraphT &G, std::string_view path) {
    std::ofstream outfile(path.data(), std::ios::binary);
    Aux::enforceOpened(outfile);
    writeData(outfile, G);
    INFO("Written graph to ", path);
}

template <class GraphT>
std::vector<uint8_t> NetworkitBinaryWriter::writeToBuffer(const GraphT &G) {
    std::stringstream outdata;
    writeData(outdata, G);
    std::string str = outdata.str();
    return std::vector<uint8_t>(str.begin(), str.end());
}

template <class StreamT, class GraphT>
void NetworkitBinaryWriter::writeData(StreamT &outStream, const GraphT &G) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    nkbg::WEIGHT_FORMAT weightFormat;

    auto detectWeightsType = [&]() -> nkbg::WEIGHT_FORMAT {
        if (!G.isWeighted())
            return nkbg::WEIGHT_FORMAT::NONE;

        if constexpr (std::is_integral_v<EdgeWeightT>) {
            return std::is_signed_v<EdgeWeightT> ? nkbg::WEIGHT_FORMAT::SIGNED_VARINT
                                                 : nkbg::WEIGHT_FORMAT::VARINT;
        }
        if constexpr (std::is_floating_point_v<EdgeWeightT>) {
            return std::is_same_v<EdgeWeightT, float> ? nkbg::WEIGHT_FORMAT::FLOAT
                                                      : nkbg::WEIGHT_FORMAT::DOUBLE;
        }
        throw std::runtime_error("Unsupported NetworkitBinaryGraph weight type");
    };

    switch (weightsType) {
    case NetworkitBinaryWeights::NONE:
        weightFormat = nkbg::WEIGHT_FORMAT::NONE;
        break;
    case NetworkitBinaryWeights::AUTO_DETECT:
        weightFormat = detectWeightsType();
        break;
    case NetworkitBinaryWeights::UNSIGNED_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::VARINT;
        break;
    case NetworkitBinaryWeights::SIGNED_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::SIGNED_VARINT;
        break;
    case NetworkitBinaryWeights::DOUBLE_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::DOUBLE;
        break;
    case NetworkitBinaryWeights::FLOAT_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::FLOAT;
        break;
    }

    switch (edgeIndex) {
    case NetworkitBinaryEdgeIDs::NO_EDGE_IDS:
        preserveEdgeIndex = false;
        break;
    case NetworkitBinaryEdgeIDs::AUTO_DETECT:
    case NetworkitBinaryEdgeIDs::WIRTE_EDGE_IDS:
        preserveEdgeIndex = G.hasEdgeIds();
        break;
    }

    nkbg::Header header;
    uint64_t adjListSize = 0;
    uint64_t adjTransposeSize = 0;

    auto setFeatures = [&](uint8_t tableWidth) {
        header.features =
            (G.isDirected() & nkbg::DIR_MASK)
            | ((static_cast<uint64_t>(weightFormat) << nkbg::WGHT_SHIFT) & nkbg::WGHT_MASK)
            | ((static_cast<uint64_t>(preserveEdgeIndex) << nkbg::INDEX_SHIFT) & nkbg::INDEX_MASK)
            | ((static_cast<uint64_t>(nkbg::widthCode(sizeof(NodeT)))
                << nkbg::NODE_TYPE_WIDTH_SHIFT)
               & nkbg::NODE_TYPE_WIDTH_MASK)
            | ((static_cast<uint64_t>(nkbg::widthCode(tableWidth)) << nkbg::TABLE_WIDTH_SHIFT)
               & nkbg::TABLE_WIDTH_MASK);
    };

    auto writeHeader = [&]() {
        outStream.write(header.magic, 8);
        outStream.write(reinterpret_cast<char *>(&header.checksum), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.features), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.nodes), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.chunks), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetBaseData), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjLists), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjTranspose), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetWeightLists), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetWeightTranspose), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjIdLists), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjIdTranspose), sizeof(uint64_t));
    };

    auto writeWeightsToFile = [&](EdgeWeightT w) {
        uint8_t tmp[10];
        uint64_t weightSize;
        switch (weightFormat) {
        case nkbg::WEIGHT_FORMAT::VARINT:
            weightSize = nkbg::varIntEncode(static_cast<uint64_t>(w), tmp);
            outStream.write(reinterpret_cast<char *>(tmp), weightSize);
            break;
        case nkbg::WEIGHT_FORMAT::DOUBLE: {
            double weight = static_cast<double>(w);
            outStream.write(reinterpret_cast<char *>(&weight), sizeof(double));
        } break;
        case nkbg::WEIGHT_FORMAT::SIGNED_VARINT: {
            const uint64_t weight = nkbg::zigzagEncode(static_cast<int64_t>(w));
            weightSize = nkbg::varIntEncode(weight, tmp);
            outStream.write(reinterpret_cast<char *>(tmp), weightSize);
        } break;
        case nkbg::WEIGHT_FORMAT::FLOAT: {
            float weight = static_cast<float>(w);
            outStream.write(reinterpret_cast<char *>(&weight), sizeof(float));
        } break;
        case nkbg::WEIGHT_FORMAT::NONE:
            break;
        }
    };

    count nodes = G.upperNodeIdBound();
    if (nodes == 0) {
        strncpy(header.magic, FILE_FORMAT, 8);
        header.checksum = 0;
        setFeatures(1);
        writeHeader();
        uint8_t zero = 0;
        outStream.write(reinterpret_cast<char *>(&zero), sizeof(zero));
        return;
    }

    if (nodes < chunks) {
        chunks = nodes;
        INFO("reducing chunks to ", chunks, " chunks");
    }

    std::vector<uint64_t> firstInChunk;
    firstInChunk.push_back(0);
    uint64_t firstNode = 0;
    for (uint64_t c = 1; c < chunks; c++) {
        firstNode += (nodes / chunks);
        firstInChunk.push_back(firstNode);
    }
    firstInChunk.push_back(nodes);

    uint64_t adjSize = 0;
    uint64_t transpSize = 0;
    uint64_t adjWeightSize = 0;
    uint64_t transpWeightSize = 0;
    uint64_t adjIndexSize = 0;
    uint64_t transpIndexSize = 0;
    std::unordered_map<NodeT, uint64_t> nrOutNbrs;
    std::unordered_map<NodeT, uint64_t> nrInNbrs;
    std::vector<uint64_t> adjOffsets;
    std::vector<uint64_t> transpOffsets;
    std::vector<uint64_t> adjWghtOffsets;
    std::vector<uint64_t> transpWghtOffsets;
    std::vector<uint64_t> adjIndexOffsets;
    std::vector<uint64_t> transpIndexOffsets;

    auto computeWeightsOffsets = [&](EdgeWeightT w) {
        uint64_t size = 0;
        uint8_t tmp[10];
        switch (weightFormat) {
        case nkbg::WEIGHT_FORMAT::VARINT:
            size = nkbg::varIntEncode(static_cast<uint64_t>(w), tmp);
            break;
        case nkbg::WEIGHT_FORMAT::DOUBLE:
            size = sizeof(double);
            break;
        case nkbg::WEIGHT_FORMAT::SIGNED_VARINT:
            size = nkbg::varIntEncode(nkbg::zigzagEncode(static_cast<int64_t>(w)), tmp);
            break;
        case nkbg::WEIGHT_FORMAT::FLOAT:
            size = sizeof(float);
            break;
        case nkbg::WEIGHT_FORMAT::NONE:
            break;
        }
        return size;
    };

    auto computeIndexOffsets = [&](edgeid id) {
        uint8_t tmp[10];
        return nkbg::varIntEncode(id, tmp);
    };

    uint64_t maxTableValue = nodes;
    for (uint64_t c = 0; c < chunks; c++) {
        for (uint64_t n = firstInChunk[c]; n < firstInChunk[c + 1]; n++) {
            uint64_t outNbrs = 0;
            uint64_t inNbrs = 0;
            uint8_t tmp[10];
            const NodeT u = static_cast<NodeT>(n);

            if (!G.hasNode(u)) {
                adjSize += nkbg::varIntEncode(0, tmp);
                transpSize += nkbg::varIntEncode(0, tmp);
                nrOutNbrs[u] = 0;
                nrInNbrs[u] = 0;
                continue;
            }
            if (!G.isDirected()) {
                G.forNeighborsOf(u, [&](NodeT v, EdgeWeightT w) {
                    if (v <= u) {
                        outNbrs++;
                        adjSize += nkbg::varIntEncode(v, tmp);
                        adjWeightSize += computeWeightsOffsets(w);
                        if (preserveEdgeIndex)
                            adjIndexSize += computeIndexOffsets(G.edgeId(u, v));
                    }
                    if (v >= u) {
                        inNbrs++;
                        transpSize += nkbg::varIntEncode(v, tmp);
                        transpWeightSize += computeWeightsOffsets(w);
                        if (preserveEdgeIndex)
                            transpIndexSize += computeIndexOffsets(G.edgeId(u, v));
                    }
                });
            } else {
                G.forNeighborsOf(u, [&](NodeT v, EdgeWeightT w) {
                    outNbrs++;
                    adjSize += nkbg::varIntEncode(v, tmp);
                    adjWeightSize += computeWeightsOffsets(w);
                    if (preserveEdgeIndex)
                        adjIndexSize += computeIndexOffsets(G.edgeId(u, v));
                });
                G.forInNeighborsOf(u, [&](NodeT v, EdgeWeightT w) {
                    inNbrs++;
                    transpSize += nkbg::varIntEncode(v, tmp);
                    transpWeightSize += computeWeightsOffsets(w);
                    if (preserveEdgeIndex)
                        transpIndexSize += computeIndexOffsets(G.edgeId(v, u));
                });
            }
            adjListSize += outNbrs;
            adjSize += nkbg::varIntEncode(outNbrs, tmp);
            nrOutNbrs[u] = outNbrs;

            adjTransposeSize += inNbrs;
            transpSize += nkbg::varIntEncode(inNbrs, tmp);
            nrInNbrs[u] = inNbrs;
        }
        adjOffsets.push_back(adjSize);
        transpOffsets.push_back(transpSize);
        adjWghtOffsets.push_back(adjWeightSize);
        transpWghtOffsets.push_back(transpWeightSize);
        adjIndexOffsets.push_back(adjIndexSize);
        transpIndexOffsets.push_back(transpIndexSize);
        maxTableValue = std::max({maxTableValue, adjSize, transpSize, adjWeightSize,
                                  transpWeightSize, adjIndexSize, transpIndexSize});
    }

    const uint8_t tableWidth = nkbg::getFitWidthBytes(maxTableValue);
    strncpy(header.magic, FILE_FORMAT, 8);
    header.checksum = 0;
    setFeatures(tableWidth);
    header.nodes = nodes;
    header.chunks = chunks;
    header.offsetBaseData = sizeof(nkbg::Header);
    auto nextOffset = header.offsetBaseData;
    nextOffset += nodes * sizeof(uint8_t) + (chunks - 1) * tableWidth;
    header.offsetAdjLists = nextOffset;
    nextOffset += (chunks - 1) * tableWidth + tableWidth + adjOffsets.back();
    header.offsetAdjTranspose = nextOffset;
    nextOffset += (chunks - 1) * tableWidth + tableWidth + transpOffsets.back();

    if (weightFormat != nkbg::WEIGHT_FORMAT::NONE) {
        header.offsetWeightLists = nextOffset;
        nextOffset += (chunks - 1) * tableWidth + adjWghtOffsets.back();
        header.offsetWeightTranspose = nextOffset;
        nextOffset += (chunks - 1) * tableWidth + transpWghtOffsets.back();
    } else {
        header.offsetWeightLists = 0;
        header.offsetWeightTranspose = 0;
        nextOffset += 2 * (chunks - 1) * tableWidth;
    }

    if (preserveEdgeIndex) {
        header.offsetAdjIdLists = nextOffset;
        nextOffset += (chunks - 1) * tableWidth + adjIndexOffsets.back();
        header.offsetAdjIdTranspose = nextOffset;
    } else {
        header.offsetAdjIdLists = 0;
        header.offsetAdjIdTranspose = 0;
    }
    writeHeader();

    for (uint64_t u = 0; u < nodes; ++u) {
        uint8_t nodeFlag = 0;
        if (!G.hasNode(static_cast<NodeT>(u)))
            nodeFlag |= nkbg::DELETED_BIT;
        outStream.write(reinterpret_cast<char *>(&nodeFlag), sizeof(uint8_t));
    }

    assert(!firstInChunk[0]);
    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, firstInChunk[c], tableWidth);

    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, adjOffsets[c - 1], tableWidth);
    nkbg::writeUint(outStream, adjListSize, tableWidth);
    for (uint64_t u64 = 0; u64 < nodes; ++u64) {
        const NodeT u = static_cast<NodeT>(u64);
        uint8_t tmp[10];
        const uint64_t nbrsSize = nkbg::varIntEncode(nrOutNbrs[u], tmp);
        outStream.write(reinterpret_cast<char *>(tmp), nbrsSize);
        if (!G.hasNode(u))
            continue;
        G.forNeighborsOf(u, [&](NodeT v) {
            if (!G.isDirected() && v > u)
                return;
            const uint64_t nodeSize = nkbg::varIntEncode(v, tmp);
            outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
        });
    }

    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, transpOffsets[c - 1], tableWidth);
    nkbg::writeUint(outStream, adjTransposeSize, tableWidth);
    for (uint64_t u64 = 0; u64 < nodes; ++u64) {
        const NodeT u = static_cast<NodeT>(u64);
        uint8_t tmp[10];
        const uint64_t nbrsSize = nkbg::varIntEncode(nrInNbrs[u], tmp);
        outStream.write(reinterpret_cast<char *>(tmp), nbrsSize);
        if (!G.hasNode(u))
            continue;
        if (!G.isDirected()) {
            G.forNeighborsOf(u, [&](NodeT v) {
                if (v < u)
                    return;
                const uint64_t nodeSize = nkbg::varIntEncode(v, tmp);
                outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
            });
        } else {
            G.forInNeighborsOf(u, [&](NodeT v) {
                const uint64_t nodeSize = nkbg::varIntEncode(v, tmp);
                outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
            });
        }
    }

    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, adjWghtOffsets[c - 1], tableWidth);
    for (uint64_t u64 = 0; u64 < nodes; ++u64) {
        const NodeT u = static_cast<NodeT>(u64);
        if (!G.hasNode(u))
            continue;
        G.forNeighborsOf(u, [&](NodeT v, EdgeWeightT w) {
            if (!G.isDirected() && v > u)
                return;
            writeWeightsToFile(w);
        });
    }

    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, transpWghtOffsets[c - 1], tableWidth);
    for (uint64_t u64 = 0; u64 < nodes; ++u64) {
        const NodeT u = static_cast<NodeT>(u64);
        if (!G.hasNode(u))
            continue;
        if (!G.isDirected()) {
            G.forNeighborsOf(u, [&](NodeT v, EdgeWeightT w) {
                if (v < u)
                    return;
                writeWeightsToFile(w);
            });
        } else {
            G.forInNeighborsOf(u, [&](NodeT, EdgeWeightT w) { writeWeightsToFile(w); });
        }
    }

    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, adjIndexOffsets[c - 1], tableWidth);
    if (preserveEdgeIndex) {
        for (uint64_t u64 = 0; u64 < nodes; ++u64) {
            const NodeT u = static_cast<NodeT>(u64);
            if (!G.hasNode(u))
                continue;
            G.forNeighborsOf(u, [&](NodeT v) {
                if (!G.isDirected() && v > u)
                    return;
                uint8_t tmp[10];
                const uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(u, v), tmp);
                outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
            });
        }
    }

    for (uint64_t c = 1; c < chunks; c++)
        nkbg::writeUint(outStream, transpIndexOffsets[c - 1], tableWidth);
    if (preserveEdgeIndex) {
        for (uint64_t u64 = 0; u64 < nodes; ++u64) {
            const NodeT u = static_cast<NodeT>(u64);
            if (!G.hasNode(u))
                continue;
            if (!G.isDirected()) {
                G.forNeighborsOf(u, [&](NodeT v) {
                    if (v < u)
                        return;
                    uint8_t tmp[10];
                    const uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(u, v), tmp);
                    outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
                });
            } else {
                G.forInNeighborsOf(u, [&](NodeT v) {
                    uint8_t tmp[10];
                    const uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(v, u), tmp);
                    outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
                });
            }
        }
    }
}

} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_WRITER_IMPL_HPP_
