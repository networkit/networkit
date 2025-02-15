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
    AttributeStorageBase(std::string name, std::type_index type)
        : name{std::move(name)}, type{type}, validStorage{true} {}

    void invalidateStorage() { validStorage = false; }

    const std::string &getName() const noexcept { return name; }

    std::type_index getType() const noexcept { return type; }

    virtual std::shared_ptr<AttributeStorageBase> clone() const = 0;
    virtual void erase(index i) = 0;
    virtual void swapData(index i, index j) = 0;

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
        if (n >= valid.size())
            valid.resize(n + 1);
        if (!valid[n]) {
            valid[n] = true;
            ++validElements;
        }
    }

    void checkIndex(index n) const {
        if (!isValid(n)) {
            throw std::runtime_error("Invalid attribute value");
        }
    }

private:
    std::string name;
    std::type_index type;

protected:
    index validElements = 0;
    bool validStorage;       // Validity of the whole storage
    std::vector<bool> valid; // For each node/edgeid: whether attribute is set or not.

}; // class AttributeStorageBase

template <typename NodeOrEdge, typename GraphType>
using ASB = AttributeStorageBase<NodeOrEdge, GraphType>;

template <typename NodeOrEdge, typename GraphType, typename T, bool isConst>
class Attribute;

template <typename NodeOrEdge, typename GraphType, template <typename, typename> class Base,
          typename T>
class AttributeStorage final : public Base<NodeOrEdge, GraphType> {
public:
    AttributeStorage(std::string name) : Base<NodeOrEdge, GraphType>{std::move(name), typeid(T)} {}

    std::shared_ptr<Base<NodeOrEdge, GraphType>> clone() const override {
        return std::make_shared<AttributeStorage>(*this);
    };

    void swapData(index i, index j) override {
        assert(i < values.size());
        assert(j < values.size());
        using std::swap;
        swap(values[i], values[j]);
        swap(this->valid[i], this->valid[j]);
    }

    void erase(index i) override {
        assert(i < values.size());
        values.erase(values.begin() + i);
        this->valid.erase(this->valid.begin() + i);
    }

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
        IndexProxy(IndexProxy &&) = delete;
        IndexProxy(const IndexProxy &) = delete;

        // reading at idx
        operator T() const {
            storage->checkIndex(idx);
            return storage->values[idx];
        }

        // writing at idx
        IndexProxy &operator=(T &&other)
            requires(!isConst)
        {
            storage->set(idx, std::move(other));
            return *this;
        }

        // move + copy assignment from other proxy objects

        IndexProxy &operator=(IndexProxy &&other) noexcept
            requires(!isConst)
        {
            storage->set(idx, static_cast<T>(other));
            return *this;
        }

        IndexProxy &operator=(const IndexProxy &other)
            requires(!isConst)
        {
            storage->set(idx, static_cast<T>(other));
            return *this;
        }

    private:
        AttributeStorage_type *storage;
        index idx;
    }; // class IndexProxy
public:
    explicit Attribute(std::shared_ptr<AttributeStorage_type> ownedStorage = nullptr,
                       const GraphType *graph = nullptr)
        : ownedStorage{ownedStorage}, theGraph{graph},
          valid{ownedStorage != nullptr && graph != nullptr} {}

    Attribute(Attribute const &other)
        : ownedStorage{other.ownedStorage}, theGraph{other.theGraph}, valid{other.valid} {}

    template <bool ic = isConst, std::enable_if_t<ic, int> = 0>
    Attribute(Attribute<NodeOrEdge, GraphType, T, false> const &other)
        : ownedStorage{other.ownedStorage}, theGraph{other.theGraph}, valid{other.valid} {}

    Attribute &operator=(Attribute other) {
        this->swap(other);
        return *this;
    }

    void swap(Attribute &other) {
        std::swap(ownedStorage, other.ownedStorage);
        std::swap(theGraph, other.theGraph);
        std::swap(valid, other.valid);
    }

    Attribute(Attribute &&other) noexcept
        : ownedStorage{std::move(other.ownedStorage)}, theGraph{std::move(other.theGraph)},
          valid{std::move(other.valid)} {
        other.valid = false;
    }

    template <bool ic = isConst, std::enable_if_t<ic, int> = 0>
    Attribute(Attribute<NodeOrEdge, GraphType, T, false> &&other) noexcept
        : ownedStorage{std::move(other.ownedStorage)}, theGraph{std::move(other.theGraph)},
          valid{std::move(other.valid)} {
        other.valid = false;
    }

    auto begin() const {
        checkAttribute();
        return Iterator(lockStorage().get()).nextValid();
    }

    auto end() const { return Iterator(nullptr); }

    auto size() const noexcept { return lockStorage()->size(); }

    template <bool ic = isConst>
    std::enable_if_t<!ic> set(index i, T v) {
        checkAttribute();
#ifndef NDEBUG
        indexOK(i);
#endif // NDEBUG
        lockStorage()->set(i, std::move(v));
    }

    template <bool ic = isConst>
    std::enable_if_t<!ic> set2(node u, node v, T t) {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        checkAttribute();
        set(theGraph->edgeId(u, v), t);
    }

    auto get(index i) const {
        checkAttribute();
#ifndef NDEBUG
        indexOK(i);
#endif // NDEBUG
        return lockStorage()->get(i);
    }

    auto get2(node u, node v) const {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        checkAttribute();
        return get(theGraph->edgeId(u, v));
    }

    auto get(index i, T defaultT) const {
        checkAttribute();
#ifndef NDEBUG
        indexOK(i);
#endif // NDEBUG
        return lockStorage()->get(i, defaultT);
    }

    auto get2(node u, node v, T defaultT) const {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        checkAttribute();
        return get(theGraph->edgeId(u, v), defaultT);
    }

    IndexProxy operator[](index i) const {
        checkAttribute();
#ifndef NDEBUG
        indexOK(i);
#endif // NDEBUG
        return IndexProxy(lockStorage().get(), i);
    }

    IndexProxy operator()(node u, node v) const {
        static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
        checkAttribute();
        auto idx = theGraph->edgeId(u, v);
#ifndef NDEBUG
        indexOK(idx);
#endif // NDEBUG
        return IndexProxy(lockStorage().get(), idx);
    }

    void checkAttribute() const {
        if (!lockStorage()->validStorage || !valid)
            throw std::runtime_error("Invalid attribute");
    }

    auto getName() const {
        checkAttribute();
        return lockStorage()->getName();
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
#ifndef NDEBUG
    void indexOK(index n) const {
        if constexpr (NodeOrEdge::edges) {
            auto uv = theGraph->edgeById(n);
            if (uv.first == none) {
                throw std::runtime_error("This edgeId does not exist");
            }
        } else {
            if (!theGraph->hasNode(n)) {
                throw std::runtime_error("This node does not exist");
            }
        }
    };
#endif // NDEBUG

    auto lockStorage() const {
        auto s = ownedStorage.lock();
        if (s)
            return s;
        throw std::runtime_error("Attribute does not exist");
    }

    std::weak_ptr<AttributeStorage_type> ownedStorage;
    const GraphType *theGraph;
    bool valid;
}; // class Attribute

template <typename NodeOrEdge, typename GraphType>
class AttributeMap {
    friend GraphType;
    const GraphType *theGraph;

    std::unordered_map<std::string, std::shared_ptr<ASB<NodeOrEdge, GraphType>>> attrMap;

public:
    AttributeMap(const GraphType *g) : theGraph{g} { assert(theGraph != nullptr); }

    // do not allow copying of AttributeMap. There is a 1:1 relation to the Graph!
    AttributeMap(AttributeMap &) = delete;
    AttributeMap &operator=(const AttributeMap &) = delete;

    // copying is only allowed with a new graph pointer. This constructor copies all data.
    AttributeMap(const AttributeMap &other, const GraphType *g) : theGraph(g) {
        assert(theGraph != nullptr);
        // manual copy is required here since the copy constructor for unordered_map would only copy
        // the shared_ptr and not the data
        for (auto &[key, value] : other.attrMap) {
            auto ptr = value->clone();
            attrMap.emplace(key, ptr);
        }
    }

    // move operators
    AttributeMap(AttributeMap &&other) = default;
    AttributeMap &operator=(AttributeMap &&other) = default;

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
        if constexpr (NodeOrEdge::edges) {
            if (!theGraph->hasEdgeIds()) {
                throw std::runtime_error("Edges must be indexed");
            }
        }
        auto ownedPtr =
            std::make_shared<AttributeStorage<NodeOrEdge, GraphType, ASB, T>>(std::string{name});
        auto insertResult = attrMap.emplace(ownedPtr->getName(), ownedPtr);
        auto success = insertResult.second;
        if (!success) {
            throw std::runtime_error("Attribute with same name already exists");
        }
        return Attribute<NodeOrEdge, GraphType, T, false>{ownedPtr, theGraph};
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
            std::static_pointer_cast<AttributeStorage<NodeOrEdge, GraphType, ASB, T>>(it->second),
            theGraph};
    }

    template <typename T>
    auto get(const std::string &name) const {
        auto it = find(name);
        if (it->second.get()->getType() != typeid(T))
            throw std::runtime_error("Type mismatch in Attributes().get()");
        return Attribute<NodeOrEdge, GraphType, T, true>{
            std::static_pointer_cast<const AttributeStorage<NodeOrEdge, GraphType, ASB, T>>(
                it->second),
            theGraph};
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
