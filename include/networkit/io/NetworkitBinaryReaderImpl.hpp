/*
 * NetworkitBinaryReaderImpl.hpp
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_READER_IMPL_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_READER_IMPL_HPP_

#include <atomic>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>

namespace NetworKit {

namespace nkbg {
inline const char *accessData(const MemoryMappedFile &source) {
    return source.cbegin();
}

inline const char *accessData(const std::vector<uint8_t> &source) {
    return reinterpret_cast<const char *>(source.data());
}
} // namespace nkbg

template <class GraphT>
GraphT NetworkitBinaryReader::read(std::string_view path) {
    MemoryMappedFile mmfile(path);
    return readData<GraphT>(mmfile);
}

template <class GraphT>
GraphT NetworkitBinaryReader::readFromBuffer(const std::vector<uint8_t> &data) {
    return readData<GraphT>(data);
}

template <class GraphT, class T>
GraphT NetworkitBinaryReader::readData(const T &source) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    nkbg::Header header;
    nkbg::WeightFormat weightFormat;

    const char *startIt = nkbg::accessData(source);
    const char *it = startIt;
    auto readHeader = [&]() {
        memcpy(&header.magic, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.checksum, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.features, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        if (!memcmp("nkbg002", header.magic, 8)) {
            version = 2;
        } else if (!memcmp("nkbg003", header.magic, 8)) {
            version = 3;
        } else if (!memcmp("nkbg004", header.magic, 8)) {
            version = 4;
        } else if (!memcmp("nkbg005", header.magic, 8)) {
            version = 5;
        } else {
            throw std::runtime_error("Reader expected another magic value");
        }
        directed = (header.features & nkbg::DIR_MASK);
        weightFormat = static_cast<nkbg::WeightFormat>((header.features & nkbg::WGHT_MASK)
                                                       >> nkbg::WGHT_SHIFT);
        indexed = false;
        if (version >= 3)
            indexed = (header.features & nkbg::INDEX_MASK) >> nkbg::INDEX_SHIFT;
        tableWidth = 8;
        if (version >= 5)
            tableWidth = nkbg::widthBytes(static_cast<uint8_t>(
                (header.features & nkbg::TABLE_WIDTH_MASK) >> nkbg::TABLE_WIDTH_SHIFT));
        memcpy(&header.nodes, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.chunks, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetBaseData, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetWeightLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetWeightTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        if (version > 2) {
            memcpy(&header.offsetAdjIdLists, it, sizeof(uint64_t));
            it += sizeof(uint64_t);
            memcpy(&header.offsetAdjIdTranspose, it, sizeof(uint64_t));
            it += sizeof(uint64_t);
        }
    };

    readHeader();
    nodes = header.nodes;
    chunks = header.chunks;
    if (nodes > static_cast<uint64_t>(std::numeric_limits<NodeT>::max()))
        throw std::runtime_error("File node id range does not fit the requested graph node type");
    weighted = weightFormat != nkbg::WeightFormat::NONE;

    GraphT G(nodes, weighted, directed);
    if (indexed)
        G.indexEdges();
    if (nodes == 0)
        return G;

    const char *baseIt = startIt + header.offsetBaseData;
    for (uint64_t i = 0; i < nodes; i++) {
        uint8_t flag;
        memcpy(&flag, baseIt, sizeof(uint8_t));
        baseIt += sizeof(uint8_t);
        if (flag & nkbg::DELETED_BIT)
            G.removeNode(static_cast<NodeT>(i));
    }

    std::vector<uint64_t> firstVert;
    firstVert.push_back(0);
    for (uint64_t ch = 1; ch < chunks; ch++) {
        firstVert.push_back(nkbg::readUint(baseIt, tableWidth));
        baseIt += tableWidth;
    }
    firstVert.push_back(nodes);

    const char *adjIt = startIt + header.offsetAdjLists;
    const char *transpIt = startIt + header.offsetAdjTranspose;
    const char *adjWghtIt = startIt + header.offsetWeightLists;
    const char *transpWghtIt = startIt + header.offsetWeightTranspose;
    const char *adjIdIt = (version > 2) ? startIt + header.offsetAdjIdLists : nullptr;
    const char *transpIdIt = (version > 2) ? startIt + header.offsetAdjIdTranspose : nullptr;
    const uint64_t adjListSize = nkbg::readUint(adjIt + (chunks - 1) * tableWidth, tableWidth);
    const uint64_t transposeListSize =
        nkbg::readUint(transpIt + (chunks - 1) * tableWidth, tableWidth);

    if (!directed)
        assert(adjListSize == transposeListSize);
    G.setEdgeCount(unsafe, adjListSize);

    std::atomic<count> selfLoops{0};
    edgeid omega = 0;

    auto constructGraph = [&](uint64_t c) {
        const NodeT vertex = static_cast<NodeT>(firstVert[c]);
        uint64_t off = 0;
        uint64_t transpOff = 0;
        uint64_t wghtOff = 0;
        uint64_t transWghtOff = 0;
        uint64_t indexOff = 0;
        uint64_t transIndexOff = 0;
        if (vertex) {
            off = nkbg::readUint(adjIt + (c - 1) * tableWidth, tableWidth);
            transpOff = nkbg::readUint(transpIt + (c - 1) * tableWidth, tableWidth);
            wghtOff = nkbg::readUint(adjWghtIt + (c - 1) * tableWidth, tableWidth);
            transWghtOff = nkbg::readUint(transpWghtIt + (c - 1) * tableWidth, tableWidth);
            if (indexed) {
                indexOff = nkbg::readUint(adjIdIt + (c - 1) * tableWidth, tableWidth);
                transIndexOff = nkbg::readUint(transpIdIt + (c - 1) * tableWidth, tableWidth);
            }
        }

        off += (chunks - 1) * tableWidth + tableWidth;
        transpOff += (chunks - 1) * tableWidth + tableWidth;
        wghtOff += (chunks - 1) * tableWidth;
        transWghtOff += (chunks - 1) * tableWidth;
        indexOff += (chunks - 1) * tableWidth;
        transIndexOff += (chunks - 1) * tableWidth;
        const uint64_t n = firstVert[c + 1] - firstVert[c];

        auto consumeAdjWeight = [&]() {
            switch (weightFormat) {
            case nkbg::WeightFormat::VARINT:
            case nkbg::WeightFormat::SIGNED_VARINT: {
                uint64_t tmp;
                wghtOff +=
                    nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjWghtIt + wghtOff), tmp);
            } break;
            case nkbg::WeightFormat::DOUBLE:
                wghtOff += sizeof(double);
                break;
            case nkbg::WeightFormat::FLOAT:
                wghtOff += sizeof(float);
                break;
            case nkbg::WeightFormat::NONE:
                break;
            }
        };

        auto consumeTranspWeight = [&]() {
            switch (weightFormat) {
            case nkbg::WeightFormat::VARINT:
            case nkbg::WeightFormat::SIGNED_VARINT: {
                uint64_t tmp;
                transWghtOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(transpWghtIt + transWghtOff), tmp);
            } break;
            case nkbg::WeightFormat::DOUBLE:
                transWghtOff += sizeof(double);
                break;
            case nkbg::WeightFormat::FLOAT:
                transWghtOff += sizeof(float);
                break;
            case nkbg::WeightFormat::NONE:
                break;
            }
        };

        auto consumeAdjId = [&]() {
            if (!indexed)
                return;
            edgeid tmpId;
            indexOff +=
                nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIdIt + indexOff), tmpId);
        };

        auto consumeTranspId = [&]() {
            if (!indexed)
                return;
            edgeid tmpId;
            transIndexOff += nkbg::varIntDecode(
                reinterpret_cast<const uint8_t *>(transpIdIt + transIndexOff), tmpId);
        };

        auto readAdjWeight = [&]() -> EdgeWeightT {
            EdgeWeightT weight{1};
            switch (weightFormat) {
            case nkbg::WeightFormat::VARINT: {
                uint64_t unsignedWeight;
                wghtOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(adjWghtIt + wghtOff), unsignedWeight);
                weight = static_cast<EdgeWeightT>(unsignedWeight);
            } break;
            case nkbg::WeightFormat::DOUBLE: {
                double doubleWeight;
                memcpy(&doubleWeight, adjWghtIt + wghtOff, sizeof(double));
                wghtOff += sizeof(double);
                weight = static_cast<EdgeWeightT>(doubleWeight);
            } break;
            case nkbg::WeightFormat::SIGNED_VARINT: {
                uint64_t unsignedWeight;
                wghtOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(adjWghtIt + wghtOff), unsignedWeight);
                weight = static_cast<EdgeWeightT>(nkbg::zigzagDecode(unsignedWeight));
            } break;
            case nkbg::WeightFormat::FLOAT: {
                float floatWeight;
                memcpy(&floatWeight, adjWghtIt + wghtOff, sizeof(float));
                wghtOff += sizeof(float);
                weight = static_cast<EdgeWeightT>(floatWeight);
            } break;
            case nkbg::WeightFormat::NONE:
                break;
            }
            return weight;
        };

        auto readTranspWeight = [&]() -> EdgeWeightT {
            EdgeWeightT weight{1};
            switch (weightFormat) {
            case nkbg::WeightFormat::VARINT: {
                uint64_t unsignedWeight;
                transWghtOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(transpWghtIt + transWghtOff), unsignedWeight);
                weight = static_cast<EdgeWeightT>(unsignedWeight);
            } break;
            case nkbg::WeightFormat::DOUBLE: {
                double doubleWeight;
                memcpy(&doubleWeight, transpWghtIt + transWghtOff, sizeof(double));
                transWghtOff += sizeof(double);
                weight = static_cast<EdgeWeightT>(doubleWeight);
            } break;
            case nkbg::WeightFormat::SIGNED_VARINT: {
                uint64_t unsignedWeight;
                transWghtOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(transpWghtIt + transWghtOff), unsignedWeight);
                weight = static_cast<EdgeWeightT>(nkbg::zigzagDecode(unsignedWeight));
            } break;
            case nkbg::WeightFormat::FLOAT: {
                float floatWeight;
                memcpy(&floatWeight, transpWghtIt + transWghtOff, sizeof(float));
                transWghtOff += sizeof(float);
                weight = static_cast<EdgeWeightT>(floatWeight);
            } break;
            case nkbg::WeightFormat::NONE:
                break;
            }
            return weight;
        };

        uint64_t ignored = 0;
        for (uint64_t i = 0; i < n; i++) {
            const NodeT curr = static_cast<NodeT>(static_cast<uint64_t>(vertex) + i);
            uint64_t outNbrs;
            off += nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIt + off), outNbrs);
            uint64_t inNbrs;
            transpOff +=
                nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(transpIt + transpOff), inNbrs);

            if (!G.hasNode(curr)) {
                for (uint64_t j = 0; j < outNbrs; ++j) {
                    off +=
                        nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIt + off), ignored);
                    consumeAdjWeight();
                    consumeAdjId();
                }
                for (uint64_t j = 0; j < inNbrs; ++j) {
                    transpOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(transpIt + transpOff), ignored);
                    consumeTranspWeight();
                    consumeTranspId();
                }
                continue;
            }

            if (!directed)
                G.preallocateUndirected(curr, outNbrs + inNbrs);
            else
                G.preallocateDirected(curr, outNbrs, inNbrs);

            for (uint64_t j = 0; j < outNbrs; j++) {
                edgeid id = 0;
                uint64_t add;
                off += nkbg::varIntDecode(reinterpret_cast<const uint8_t *>(adjIt + off), add);
                const EdgeWeightT weight = readAdjWeight();
                if (indexed) {
                    indexOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(adjIdIt + indexOff), id);
                    if (id > omega)
                        omega = id;
                }
                const NodeT addNode = static_cast<NodeT>(add);
                if (!directed) {
                    if (!G.addPartialEdge(unsafe, curr, addNode, weight, id, true))
                        WARN("Not adding edge ", curr, "-", addNode,
                             " since it is already present.");
                } else {
                    if (!G.addPartialOutEdge(unsafe, curr, addNode, weight, id, true))
                        WARN("Not adding edge ", curr, "-", addNode,
                             " since it is already present.");
                }
                if (curr == addNode)
                    selfLoops.fetch_add(1, std::memory_order_relaxed);
            }

            for (uint64_t j = 0; j < inNbrs; j++) {
                uint64_t add;
                edgeid id = 0;
                transpOff += nkbg::varIntDecode(
                    reinterpret_cast<const uint8_t *>(transpIt + transpOff), add);
                const EdgeWeightT weight = readTranspWeight();
                if (indexed) {
                    transIndexOff += nkbg::varIntDecode(
                        reinterpret_cast<const uint8_t *>(transpIdIt + transIndexOff), id);
                    if (id > omega)
                        omega = id;
                }
                const NodeT addNode = static_cast<NodeT>(add);
                if (!directed) {
                    if (curr != addNode) {
                        if (G.addPartialEdge(unsafe, curr, addNode, weight, id, true))
                            WARN("Not adding edge ", curr, "-", addNode,
                                 " since it is already present.");
                    }
                } else {
                    if (!G.addPartialInEdge(unsafe, curr, addNode, weight, id, true))
                        WARN("Not adding edge ", curr, "-", addNode,
                             " since it is already present.");
                }
            }
        }
    };

#pragma omp parallel for
    for (omp_index c = 0; c < static_cast<omp_index>(chunks); c++)
        constructGraph(c);

    G.setNumberOfSelfLoops(unsafe, selfLoops.load(std::memory_order_relaxed));
    if (indexed)
        G.setUpperEdgeIdBound(unsafe, omega);
    return G;
}

} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_READER_IMPL_HPP_
