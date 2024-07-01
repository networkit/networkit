/*
 * Attributes.hpp
 *
 *  Created on: 14.05.2024
 *      Author: Klaus Ahrens
 *              Eugenio Angriman
 *              Lukas Berner
 *              Fabian Brandt-Tumescheit
 *              Alexander van der Grinten
 */

#ifndef NETWORKIT_GRAPH_ATTRIBUTES_HPP_
#define NETWORKIT_GRAPH_ATTRIBUTES_HPP_

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

// base class for all node (and edge) attribute
// storages with attribute type info
// independent of the attribute type, holds bookkeeping info only:
// - attribute name
// - type info of derived (real storage holding) classes
// - which indices are valid
// - number of valid indices
// - the associated graph (who knows, which nodes/edges exist)
// - the validity of the whole storage (initially true, false after detach)
// all indexed accesses by NetworKit::index: synonym both for node and edgeid

template <typename NodeOrEdge, typename GraphType>
class AttributeStorageBase { // alias ASB
public:
    AttributeStorageBase(const GraphType *graph, std::string name, std::type_index type)
        : name{std::move(name)}, type{type}, theGraph{graph}, validStorage{true} {
        checkPremise(); // node for PerNode, theGraph.hasEdgeIds() for PerEdges
    }

    void invalidateStorage() { validStorage = false; }

    const std::string &getName() const noexcept { return name; }

    std::type_index getType() const noexcept { return type; }

    bool isValid(index n) const noexcept { return n < valid.size() && valid[n]; }

    // Called by Graph when node/edgeid n is deleted.
    void invalidate(index n) {
        if (isValid(n)) {
            valid[n] = false;
            --validElements;
        }
    }

protected:
    void markValid(index n) {
        indexOK(n); // specialized for node/edgeid
        if (n >= valid.size())
            valid.resize(n + 1);
        if (!valid[n]) {
            valid[n] = true;
            ++validElements;
        }
    }

    void checkIndex(index n) const {
        indexOK(n);
        if (!isValid(n)) {
            throw std::runtime_error("Invalid attribute value");
        }
    }

private:
    std::string name;
    std::type_index type;
    std::vector<bool> valid; // For each node/edgeid: whether attribute is set or not.

protected:
    void indexOK(index n) const;
    void checkPremise() const;
    index validElements = 0;
    const GraphType *theGraph;
    bool validStorage; // Validity of the whole storage

}; // class AttributeStorageBase

template <typename NodeOrEdge, typename GraphType>
using ASB = AttributeStorageBase<NodeOrEdge, GraphType>;

template <typename NodeOrEdge, typename GraphType, typename T, bool isConst>
class Attribute;

template <typename NodeOrEdge, typename GraphType, template <typename, typename> class Base,
          typename T>
class AttributeStorage : public Base<NodeOrEdge, GraphType> {
public:
    AttributeStorage(const GraphType *theGraph, std::string name)
        : Base<NodeOrEdge, GraphType>{theGraph, std::move(name), typeid(T)} {}

    void resize(index i) {
        if (i >= values.size())
            values.resize(i + 1);
    }

    auto size() const noexcept { return this->validElements; }

    void set(index i, T &&v) {
        this->markValid(i);
        resize(i);
        values[i] = std::move(v);
    }

    // instead of returning an std::optional (C++17) we provide these
    // C++14 options
    // (1) throw an exception when invalid:
    T get(index i) const { // may throw
        this->checkIndex(i);
        return values[i];
    }

    // (2) give default value when invalid:
    T get(index i, T defaultT) const noexcept {
        if (i >= values.size() || !this->isValid(i))
            return defaultT;
        return values[i];
    }

    friend Attribute<NodeOrEdge, GraphType, T, true>;
    friend Attribute<NodeOrEdge, GraphType, T, false>;

private:
    using Base<NodeOrEdge, GraphType>::theGraph;
    std::vector<T> values; // the real attribute storage
};                         // class AttributeStorage<NodeOrEdge, Base, T>

template <typename NodeOrEdge, typename GraphType, typename T, bool isConst>
class Attribute {
public:
    using AttributeStorage_type =
        std::conditional_t<isConst, const AttributeStorage<NodeOrEdge, GraphType, ASB, T>,
                           AttributeStorage<NodeOrEdge, GraphType, ASB, T>>;
    class Iterator {
    public:
        // The value type of the attribute. Returned by
        // operator*().
        using value_type = T;

        // Reference to the value_type, required by STL.
        using reference = std::conditional_t<isConst, const value_type &, value_type &>;

        // Pointer to the value_type, required by STL.
        using pointer = std::conditional_t<isConst, const value_type *, value_type *>;

        // STL iterator category.
        using iterator_category = std::forward_iterator_tag;

        // Signed integer type of the result of subtracting two pointers,
        // required by STL.
        using difference_type = ptrdiff_t;

        Iterator() : storage{nullptr}, idx{0} {}
        Iterator(AttributeStorage_type *storage) : storage{storage}, idx{0} {
            if (storage) {
                nextValid();
            }
        }

        Iterator &nextValid() {
            while (storage && !storage->isValid(idx)) {
                if (idx >= storage->values.size()) {
                    storage = nullptr;
                    return *this;
                }
                ++idx;
            }
            return *this;
        }

        Iterator &operator++() {
            if (!storage) {
                throw std::runtime_error("Invalid attribute iterator");
            }
            ++idx;
            return nextValid();
        }

        auto operator*() const {
            if (!storage) {
                throw std::runtime_error("Invalid attribute iterator");
            }
            return std::make_pair(idx, storage->values[idx]);
        }

        bool operator==(Iterator const &iter) const noexcept {
            if (storage == nullptr && iter.storage == nullptr) {
                return true;
            }
            return storage == iter.storage && idx == iter.idx;
        }

        bool operator!=(Iterator const &iter) const noexcept { return !(*this == iter); }

    private:
        AttributeStorage_type *storage;
        index idx;
    }; // class Iterator

private:
    class IndexProxy {
        // a helper class for distinguished read and write on an indexed
        // attribute
        // operator[] on an attribute yields an IndexProxy holding
        // location and index of access
        //    - casting an IndexProxy to the attribute type reads the value
        //    - assigning to it (operator=) writes the value
    public:
        IndexProxy(AttributeStorage_type *storage, index idx) : storage{storage}, idx{idx} {}

        // reading at idx
        operator T() const {
            storage->checkIndex(idx);
            return storage->values[idx];
        }

        // writing at idx
        template <bool ic = isConst>
        std::enable_if_t<!ic, T> &operator=(T &&other) {
            storage->set(idx, std::move(other));
            return storage->values[idx];
        }

    private:
        AttributeStorage_type *storage;
        index idx;
    }; // class IndexProxy
public:
    explicit Attribute(std::shared_ptr<AttributeStorage_type> ownedStorage = nullptr)
        : ownedStorage{ownedStorage}, valid{ownedStorage != nullptr} {}

    Attribute(Attribute const &other) : ownedStorage{other.ownedStorage}, valid{other.valid} {}

    template <bool ic = isConst, std::enable_if_t<ic, int> = 0>
    Attribute(Attribute<NodeOrEdge, GraphType, T, false> const &other)
        : ownedStorage{other.ownedStorage}, valid{other.valid} {}

    Attribute &operator=(Attribute other) {
        this->swap(other);
        return *this;
    }

    void swap(Attribute &other) {
        std::swap(ownedStorage, other.ownedStorage);
        std::swap(valid, other.valid);
    }

    Attribute(Attribute &&other) noexcept
        : ownedStorage{std::move(other.ownedStorage)}, valid{other.valid} {
        other.valid = false;
    }

    template <bool ic = isConst, std::enable_if_t<ic, int> = 0>
    Attribute(Attribute<NodeOrEdge, GraphType, T, false> &&other) noexcept
        : ownedStorage{std::move(other.ownedStorage)}, valid{other.valid} {
        other.valid = false;
    }

    auto begin() const {
        checkAttribute();
        return Iterator(ownedStorage.get()).nextValid();
    }

    auto end() const { return Iterator(nullptr); }

    auto size() const noexcept { return ownedStorage->size(); }

    template <bool ic = isConst>
    std::enable_if_t<!ic> set(index i, T v) {
        checkAttribute();
        ownedStorage->set(i, std::move(v));
    }

    template <bool ic = isConst>
    std::enable_if_t<!ic> set2(node u, node v, T t) {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        set(ownedStorage->theGraph->edgeId(u, v), t);
    }

    auto get(index i) const {
        checkAttribute();
        return ownedStorage->get(i);
    }

    auto get2(node u, node v) const {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        return get(ownedStorage->theGraph->edgeId(u, v));
    }

    auto get(index i, T defaultT) const {
        checkAttribute();
        return ownedStorage->get(i, defaultT);
    }

    auto get2(node u, node v, T defaultT) const {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        return get(ownedStorage->theGraph->edgeId(u, v), defaultT);
    }

    IndexProxy operator[](index i) const {
        checkAttribute();
        return IndexProxy(ownedStorage.get(), i);
    }

    IndexProxy operator()(node u, node v) const {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        checkAttribute();
        return IndexProxy(ownedStorage.get(), ownedStorage->theGraph->edgeId(u, v));
    }

    void checkAttribute() const {
        if (!ownedStorage->validStorage)
            throw std::runtime_error("Invalid attribute");
    }

    auto getName() const {
        checkAttribute();
        return ownedStorage->getName();
    }

    void write(std::string const &filename) const {
        std::ofstream out(filename);
        if (!out)
            ERROR("cannot open ", filename, " for writing");

        for (auto it = begin(); it != end(); ++it) {
            auto pair = *it;
            auto n = pair.first;  // node/edgeid
            auto v = pair.second; // value
            out << n << "\t" << v << "\n";
        }
        out.close();
    }

    template <bool ic = isConst>
    std::enable_if_t<!ic> read(const std::string &filename) {
        std::ifstream in(filename);
        if (!in) {
            ERROR("cannot open ", filename, " for reading");
        }
        index n; // node/edgeid
        T v;     // value
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream istring(line);
            if constexpr (std::is_same_v<T, std::string>) {
                istring >> n >> std::ws;
                std::getline(istring, v);
            } else {
                istring >> n >> v;
            }
            set(n, v);
        }
    }

private:
    std::shared_ptr<AttributeStorage_type> ownedStorage;
    bool valid;
}; // class Attribute

template <typename NodeOrEdge, typename GraphType>
class AttributeMap {
    friend GraphType;
    const GraphType *theGraph;

public:
    std::unordered_map<std::string, std::shared_ptr<ASB<NodeOrEdge, GraphType>>> attrMap;

    AttributeMap(const GraphType *g) : theGraph{g} {}

    auto find(std::string const &name) {
        auto it = attrMap.find(name);
        if (it == attrMap.end()) {
            throw std::runtime_error("No such attribute");
        }
        return it;
    }

    auto find(std::string const &name) const {
        auto it = attrMap.find(name);
        if (it == attrMap.end()) {
            throw std::runtime_error("No such attribute");
        }
        return it;
    }

    template <typename T>
    auto attach(const std::string &name) {
        auto ownedPtr = std::make_shared<AttributeStorage<NodeOrEdge, GraphType, ASB, T>>(
            theGraph, std::string{name});
        auto insertResult = attrMap.emplace(ownedPtr->getName(), ownedPtr);
        auto success = insertResult.second;
        if (!success) {
            throw std::runtime_error("Attribute with same name already exists");
        }
        return Attribute<NodeOrEdge, GraphType, T, false>{ownedPtr};
    }

    void detach(const std::string &name) {
        auto it = find(name);
        auto storage = it->second.get();
        storage->invalidateStorage();
        it->second.reset();
        attrMap.erase(name);
    }

    template <typename T>
    auto get(const std::string &name) {
        auto it = find(name);
        if (it->second.get()->getType() != typeid(T))
            throw std::runtime_error("Type mismatch in Attributes().get()");
        return Attribute<NodeOrEdge, GraphType, T, false>{
            std::static_pointer_cast<AttributeStorage<NodeOrEdge, GraphType, ASB, T>>(it->second)};
    }

    template <typename T>
    auto get(const std::string &name) const {
        auto it = find(name);
        if (it->second.get()->getType() != typeid(T))
            throw std::runtime_error("Type mismatch in Attributes().get()");
        return Attribute<NodeOrEdge, GraphType, T, true>{
            std::static_pointer_cast<const AttributeStorage<NodeOrEdge, GraphType, ASB, T>>(
                it->second)};
    }

}; // class AttributeMap

/// @private
class PerNode {
public:
    static constexpr bool edges = false;
};

/// @private
class PerEdge {
public:
    static constexpr bool edges = true;
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_ATTRIBUTES_HPP_
