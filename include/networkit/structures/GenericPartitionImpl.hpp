#ifndef NETWORKIT_STRUCTURES_GENERIC_PARTITION_IMPL_HPP_
#define NETWORKIT_STRUCTURES_GENERIC_PARTITION_IMPL_HPP_

#include <algorithm>
#include <atomic>
#include <memory>

#include <networkit/auxiliary/Parallel.hpp>
namespace NetworKit {

template <IntegralValue IndexType>
GenericPartition<IndexType>::GenericPartition() : z(0), omega(0), data(0) {}

template <IntegralValue IndexType>
GenericPartition<IndexType>::GenericPartition(const std::vector<IndexType> &data)
    : z(data.size()), omega(), data(data) {
    auto max_elem = *std::max_element(data.begin(), data.end());
    this->omega = (max_elem == noneIndex) ? 0 : max_elem;
}

template <IntegralValue IndexType>
GenericPartition<IndexType>::GenericPartition(IndexType z) : z(z), omega(0), data(z, noneIndex) {}

template <IntegralValue IndexType>
GenericPartition<IndexType>::GenericPartition(IndexType z, IndexType defaultValue)
    : z(z), omega(0), data(z, defaultValue) {}

template <IntegralValue IndexType>
void GenericPartition<IndexType>::allToSingletons() {
    setUpperBound(numberOfElements());
    parallelForEntries([&](IndexType e, IndexType) { data[e] = e; });
}

template <IntegralValue IndexType>
IndexType GenericPartition<IndexType>::mergeSubsets(IndexType s, IndexType t) {
    assert(s <= omega);
    assert(t <= omega);
    if (s != t) {
        IndexType m = newSubsetId(); // new id for merged set
        for (IndexType e = 0; e < this->z; ++e) {
            if (data[e] == s || data[e] == t) {
                data[e] = m;
            }
        }
        return m;
    }
    return noneIndex; // no new cluster formed
}

template <IntegralValue IndexType>
count GenericPartition<IndexType>::numberOfSubsets() const {
    auto n = upperBound();
    // a boolean vector would not be thread-safe
    std::unique_ptr<std::atomic<bool>[]> exists(new std::atomic<bool>[n] {});
    this->parallelForEntries([&](IndexType, IndexType s) {
        if (s != noneIndex) {
            exists[s] = true;
        }
    });
    count k = 0; // number of actually existing clusters
#pragma omp parallel for reduction(+ : k)
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        if (exists[i]) {
            k++;
        }
    }
    return k;
}

template <IntegralValue IndexType>
void GenericPartition<IndexType>::compact(bool useTurbo) {
    IndexType i = 0;
    if (!useTurbo) {
        std::vector<IndexType> usedIds(data);
        Aux::Parallel::sort(usedIds.begin(), usedIds.end());
        auto last = std::unique(usedIds.begin(), usedIds.end());
        usedIds.erase(last, usedIds.end());
        usedIds.erase(std::remove(usedIds.begin(), usedIds.end(), noneIndex), usedIds.end());
        i = usedIds.size();

        this->parallelForEntries(
            [&](IndexType e, IndexType s) { // replace old SubsetIDs with the new IDs
                if (s != noneIndex) {
                    data[e] = std::distance(usedIds.begin(),
                                            std::lower_bound(usedIds.begin(), usedIds.end(), s));
                }
            });
    } else {
        std::vector<IndexType> compactingMap(this->upperBound(), noneIndex);
        this->forEntries([&](IndexType, IndexType s) {
            if (s != noneIndex && compactingMap[s] == noneIndex) {
                compactingMap[s] = i++;
            }
        });
        this->parallelForEntries(
            [&](IndexType e, IndexType s) { // replace old SubsetIDs with the new IDs
                if (s != noneIndex) {
                    data[e] = compactingMap[s];
                }
            });
    }
    this->setUpperBound(i);
}

template <IntegralValue IndexType>
std::vector<count> GenericPartition<IndexType>::subsetSizes() const {
    std::vector<count> sizes;
    std::map<IndexType, count> map = this->subsetSizeMap();
    sizes.reserve(map.size());
    for (auto kv : map) {
        sizes.push_back(kv.second);
    }
    return sizes;
}

template <IntegralValue IndexType>
std::map<IndexType, count> GenericPartition<IndexType>::subsetSizeMap() const {
    std::map<IndexType, count> subset2size;

    this->forEntries([&](IndexType, IndexType s) {
        if (s != noneIndex) {
            subset2size[s] += 1;
        }
    });

    return subset2size;
}

template <IntegralValue IndexType>
std::set<IndexType> GenericPartition<IndexType>::getMembers(IndexType s) const {
    assert(s <= omega);
    std::set<IndexType> subset;
    for (IndexType e = 0; e < this->z; ++e) {
        if (data[e] == s) {
            subset.insert(e);
        }
    }
    return subset;
}

template <IntegralValue IndexType>
const std::vector<IndexType> &GenericPartition<IndexType>::getVector() const {
    return data;
}

template <IntegralValue IndexType>
std::set<std::set<IndexType>> GenericPartition<IndexType>::getSubsets() const {
    std::vector<std::set<IndexType>> table(omega + 1);
    this->forEntries([&](IndexType e, IndexType s) {
        assert(s <= omega);
        table[s].insert(e);
    });

    std::set<std::set<IndexType>> subsets;
    for (const auto &set : table) {
        if (!set.empty()) {
            subsets.insert(set);
        }
    }
    return subsets;
}

template <IntegralValue IndexType>
void GenericPartition<IndexType>::allToOnePartition() {
    omega = 0;
    this->parallelForEntries([&](IndexType e, IndexType) { this->data[e] = 0; });
}

template <IntegralValue IndexType>
std::set<IndexType> GenericPartition<IndexType>::getSubsetIds() const {
    std::set<IndexType> ids;
    for (IndexType id : data) {
        if (id != noneIndex) {
            ids.insert(id);
        }
    }
    return ids;
}

} // namespace NetworKit
#endif // NETWORKIT_STRUCTURES_GENERIC_PARTITION_IMPL_HPP_
