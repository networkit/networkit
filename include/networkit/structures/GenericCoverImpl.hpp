#ifndef NETWORKIT_STRUCTURES_GENERIC_COVER_IMPL_HPP_
#define NETWORKIT_STRUCTURES_GENERIC_COVER_IMPL_HPP_

#include <algorithm>
#include <cassert>
#include <iterator>
#include <unordered_set>

namespace NetworKit {

template <IntegralValue IndexType>
GenericCover<IndexType>::GenericCover() : z(0), omega(0), data(0) {}

template <IntegralValue IndexType>
GenericCover<IndexType>::GenericCover(IndexType z) : z(z - 1), omega(0), data(z) {}

template <IntegralValue IndexType>
GenericCover<IndexType>::GenericCover(const GenericPartition<IndexType> &p)
    : z(p.numberOfElements() - 1), omega(p.upperBound() - 1), data(p.numberOfElements()) {
    p.forEntries([&](IndexType e, IndexType s) {
        if (s != GenericPartition<IndexType>::noneIndex)
            data[e].insert(s);
    });
}

template <IntegralValue IndexType>
bool GenericCover<IndexType>::contains(IndexType e) const {
    // e is in the element index range and the entry is not empty
    return (e <= z) && (!data[e].empty());
}

template <IntegralValue IndexType>
bool GenericCover<IndexType>::inSameSubset(IndexType e1, IndexType e2) const {
    assert(e1 <= z);
    assert(e2 <= z);
    assert(!data[e1].empty());
    // elements cannot be unassigned - it may be possible to change this behavior
    assert(!data[e2].empty());
    std::unordered_set<IndexType> intersect;
    std::set_intersection(data[e1].begin(), data[e1].end(), data[e2].begin(), data[e2].end(),
                          std::inserter(intersect, intersect.begin()));
    return (!intersect.empty());
}

template <IntegralValue IndexType>
std::set<IndexType> GenericCover<IndexType>::getMembers(IndexType s) const {
    assert(s <= omega);
    std::set<IndexType> members;
    for (IndexType e = 0; e <= this->z; ++e) {
        for (IndexType t : data[e]) {
            if (t == s) {
                members.insert(e);
            }
        }
    }
    return members;
}

template <IntegralValue IndexType>
void GenericCover<IndexType>::addToSubset(IndexType s, IndexType e) {
    assert(e <= z);
    assert(s <= omega);
    data[e].insert(s);
}

template <IntegralValue IndexType>
void GenericCover<IndexType>::removeFromSubset(IndexType s, IndexType e) {
    assert(e <= z);
    assert(s <= omega);
    data[e].erase(s);
}

template <IntegralValue IndexType>
void GenericCover<IndexType>::moveToSubset(IndexType s, IndexType e) {
    assert(e <= z);
    assert(s <= omega);
    data[e].clear();
    data[e].insert(s);
}

template <IntegralValue IndexType>
IndexType GenericCover<IndexType>::toSingleton(IndexType e) {
    assert(e <= z);
    data[e].clear();
    IndexType sid = newSubsetId();
    data[e].insert(sid);
    return sid;
}

template <IntegralValue IndexType>
void GenericCover<IndexType>::allToSingletons() {
    for (IndexType e = 0; e <= this->z; ++e) {
        toSingleton(e);
    }
}

template <IntegralValue IndexType>
void GenericCover<IndexType>::mergeSubsets(IndexType s, IndexType t) {
    assert(s <= omega);
    assert(t <= omega);
    if (s != t) {
        IndexType m = newSubsetId(); // new id for merged set
        for (IndexType e = 0; e <= this->z; ++e) {
            auto its = data[e].find(s);
            auto itt = data[e].find(t);
            if (its != data[e].end()) {
                data[e].erase(its);
                data[e].insert(m);
            }
            // was else if. makes errors, in case an element is in s as well as t
            if (itt != data[e].end()) {
                data[e].erase(itt);
                data[e].insert(m);
            }
        }
    }
}

template <IntegralValue IndexType>
IndexType GenericCover<IndexType>::upperBound() const {
    return omega + IndexType{1}; // to enable usual loop test x < upperBound()
}

template <IntegralValue IndexType>
IndexType GenericCover<IndexType>::lowerBound() const {
    return IndexType{0};
}

template <IntegralValue IndexType>
std::vector<count> GenericCover<IndexType>::subsetSizes() const {
    std::map<IndexType, count> mapping;
    std::vector<count> sizes;
    count newIndex = 0;
    for (IndexType e = 0; e <= this->z; ++e) { // stores sizes in a vector
        for (IndexType t : data[e]) {
            if (mapping.find(t) == mapping.end()) {
                mapping[t] = newIndex++;
                sizes.push_back(1);
            } else {
                sizes[mapping[t]]++;
            }
        }
    }
    return sizes;
}

template <IntegralValue IndexType>
std::map<IndexType, count> GenericCover<IndexType>::subsetSizeMap() const {
    std::map<IndexType, count> sizeMap;
    for (IndexType e = 0; e <= this->z; ++e) { // stores sizes of subsets in a map
        for (IndexType t : data[e]) {
            if (sizeMap.find(t) == sizeMap.end()) {
                sizeMap[t] = 1;
            } else {
                sizeMap[t]++;
            }
        }
    }
    return sizeMap;
}

template <IntegralValue IndexType>
count GenericCover<IndexType>::numberOfSubsets() const {
    std::vector<int> exists(upperBound(), 0); // a boolean vector would not be thread-safe

    this->parallelForEntries([&](IndexType, const std::set<IndexType> &s) {
        if (!s.empty()) {
            for (auto it = s.begin(); it != s.end(); it++) {
                IndexType currentSubset = *it;
                exists[currentSubset] = 1;
            }
        }
    });

    count k = 0; // number of actually existing clusters
#pragma omp parallel for reduction(+ : k)
    for (omp_index i = 0; i < static_cast<omp_index>(upperBound()); ++i) {
        if (exists[i]) {
            k++;
        }
    }

    return k;
}

template <IntegralValue IndexType>
count GenericCover<IndexType>::numberOfElements() const {
    return static_cast<count>(z + 1);
}

template <IntegralValue IndexType>
IndexType GenericCover<IndexType>::extend() {
    data.emplace_back();
    ++z;
    assert(z + 1 == data.size());
    return z;
}

template <IntegralValue IndexType>
void GenericCover<IndexType>::setUpperBound(IndexType upper) {
    this->omega = upper - 1;
}

template <IntegralValue IndexType>
std::set<IndexType> GenericCover<IndexType>::getSubsetIds() const {
    std::set<IndexType> ids;
    for (const auto &subset : data) {
        ids.insert(subset.begin(), subset.end());
    }
    return ids;
}

} /* namespace NetworKit */

#endif // NETWORKIT_STRUCTURES_GENERIC_COVER_IMPL_HPP_
