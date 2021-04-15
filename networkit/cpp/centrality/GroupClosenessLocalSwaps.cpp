// networkit-format

#include <algorithm>
#include <cassert>
#include <omp.h>
#include <queue>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/GroupClosenessLocalSwaps.hpp>

namespace NetworKit {

GroupClosenessLocalSwaps::GroupClosenessLocalSwaps(const Graph &G, const std::vector<node> &group,
                                                   count maxSwaps)
    : GroupClosenessLocalSwaps(G, group.begin(), group.end(), maxSwaps) {}

void GroupClosenessLocalSwaps::init() {
    const auto n = G->upperNodeIdBound();

    distance.assign(n, 0);
    gamma.assign(n * group.size(), false);
    canSwap.assign(group.size(), false);
    visited.assign(n, false);

    idxMap.clear();
    idxMap.reserve(group.size());

    stack.assign(n, 0);
    value.assign(group.size(), 0);
    valueDecrement.assign(group.size(), 0);
    sumOfMins.assign(n, 0);

    intDistributions.resize(omp_get_max_threads());
    urngs.reserve(omp_get_max_threads());
    for (index i = 0; i < static_cast<index>(omp_get_max_threads()); ++i)
        urngs.emplace_back(Aux::Random::getURNG());

    for (size_t i = 0; i < group.size(); ++i) {
        const auto u = group[i];
        idxMap[u] = i;
        gamma[u * group.size() + i] = 1;
    }

#ifdef AVX2_AVAILABLE
    randVec.resize(n);
#else
    randVec.resize(K * n);
#endif // AVX2_AVAILABLE

    totalSwaps = 0;
    hasRun = false;
}

void GroupClosenessLocalSwaps::run() {
    init();

    while (findAndSwap() && ++totalSwaps < maxSwaps) {
    }

    hasRun = true;
}

bool GroupClosenessLocalSwaps::findAndSwap() {
    bfsFromGroup();
    const auto v = estimateHighestDecrement();
    const auto farnessDecrement = computeFarnessDecrement(v);
    int64_t improvement = 0;
    node u = none;

    G->forNeighborsOf(v, [&](const node y) {
        if (distance[y] == 0) {
            const auto idx = idxMap.at(y);
            const auto curImprovement = farnessDecrement - value[idx] + valueDecrement[idx];

            if (curImprovement > improvement) {
                improvement = curImprovement;
                u = y;
            }
        }
    });

    for (const auto &idx : idxMap) {
        if (canSwap[idx.second]) {
            int64_t curImprovement =
                farnessDecrement - value[idx.second] + valueDecrement[idx.second];
            if (curImprovement > improvement) {
                improvement = curImprovement;
                u = idx.first;
            }
        }
    }

    if (improvement <= 0)
        return false;

    const auto idxU = idxMap.at(u);
    idxMap.erase(u);
    idxMap[v] = idxU;
    resetGamma(v, idxU);

    return true;
}

void GroupClosenessLocalSwaps::bfsFromGroup() {
    std::queue<node> q;
    std::fill(visited.begin(), visited.end(), false);
    stackSize = 0;

    for (const auto &idx : idxMap) {
        q.push(idx.first);
        value[idx.second] = 1;
        distance[idx.first] = 0;
        visited[idx.first] = true;
    }

    do {
        const auto u = q.front();
        q.pop();

        // Cannot swap with nodes in the group
        if (distance[u] == 0)
            canSwap[idxMap.at(u)] = 0;

        bool uIsLeaf = false;

        G->forNeighborsOf(u, [&](const node v) {
            // Whether v is in \Gamma_u i.e., the shortest path from S to v is realized only
            // by u.
            bool inGamma = true;

            // Whether the node in the group that realizes the shortest distance to v has
            // been found.
            bool nearestNodeFound = false;

            // Index of the node in the group that realizes the shortest distance to v.
            index groupIdx;

            if (!visited[v]) {
                uIsLeaf = false;
                distance[v] = distance[u] + 1;
                visited[v] = true;
                q.push(v);

                for (size_t i = 0; i < group.size(); ++i) {
                    const auto curGamma = gamma[group.size() * u + i];
                    if (curGamma) {
                        if (!nearestNodeFound) {
                            nearestNodeFound = true;
                            groupIdx = i;
                        } else
                            inGamma = false;
                    }

                    gamma[group.size() * v + i] = curGamma;
                }

                if (inGamma)
                    ++value[groupIdx];

            } else if (distance[u] + 1 == distance[v]) {
                inGamma = true;
                nearestNodeFound = false;
                bool subtract = false;

                for (size_t i = 0; i < group.size(); ++i) {
                    if (gamma[group.size() * v + i]) {
                        if (!nearestNodeFound) {
                            nearestNodeFound = true;
                            groupIdx = i;
                        } else {
                            inGamma = false;
                            break;
                        }
                    } else if (gamma[group.size() * u + i]) {
                        gamma[group.size() * v + i] = 1;
                        subtract = true;
                    }
                }
                if (inGamma && subtract)
                    --value[groupIdx];
            }

            if (!distance[u] && !distance[v])
                canSwap[idxMap[u]] = true;
        });

        if (distance[u] != 0 && (!uIsLeaf || distance[u] == count{1}))
            stack[stackSize++] = u;

    } while (!q.empty());
}

int64_t GroupClosenessLocalSwaps::computeFarnessDecrement(node v) {
    std::fill(visited.begin(), visited.end(), false);

    std::queue<node> q;
    q.push(v);
    distance[v] = 0;
    visited[v] = true;

    int16_t decrement{1};
    std::fill(valueDecrement.begin(), valueDecrement.end(), int64_t{0});

    do {
        const auto u = q.front();
        q.pop();

        bool inGamma = false;
        index groupIdx;

        for (size_t i = 0; i < group.size(); ++i) {
            if (gamma[group.size() * u + i]) {
                if (!inGamma) {
                    inGamma = true;
                    groupIdx = i;
                } else {
                    inGamma = false;
                    break;
                }
            }
        }

        if (inGamma)
            ++valueDecrement[groupIdx];

        G->forNeighborsOf(u, [&](const node v) {
            if (visited[v])
                return;

            if (distance[u] + 1 <= distance[v]) {
                if (distance[u] + 1 < distance[v]) {
                    distance[v] = distance[u] + 1;
                    ++decrement;
                }
                q.push(v);
            }
            visited[v] = true;
        });

    } while (!q.empty());

    return decrement;
}

void GroupClosenessLocalSwaps::initRandomVector() {
    G->parallelForNodes([&](const node u) {
        // Avoid to generate numbers for nodes in the group
        if (distance[u] > 0) {
            auto tid = omp_get_thread_num();
            auto &curUrng = urngs[tid].get();
            auto &distr = intDistributions[tid];
#ifdef AVX2_AVAILABLE
            // Generating two 16-bit random integers per time
            for (index j = 0; j < K; j += 2) {
                const auto x = distr(curUrng);
                randVec[u].items[j] = static_cast<uint16_t>(x);
                randVec[u].items[j + 1] = static_cast<uint16_t>(x >> K);
            }

            randVec[u].vec = *(__m256i *)(&randVec[u].items[0]);
#else
            // Generating two 16-bit random integers per time
            for (index j = 0; j < K; j += 2) {
                const auto x = distr(curUrng);
                randVec[K * u + j] = static_cast<uint16_t>(x);
                randVec[K * u + j + 1] = static_cast<uint16_t>(x >> K);
            }
#endif // AVX2_AVAILABLE
        }
    });
}

node GroupClosenessLocalSwaps::estimateHighestDecrement() {
    initRandomVector();

    float bestEstimate = -1.f;
    node v = none;

    for (count i = 0; i < stackSize; ++i) {
        const auto x = stack[stackSize - 1 - i];
#ifdef AVX2_AVAILABLE
        // 16 randomly generated integers;
        __m256i &x1 = randVec[x].vec;
        // Pulling leaves
        G->forNeighborsOf(x, [&](const node y) {
            if (distance[y] == distance[x] + 1) {
                const __m256i &y1 = randVec[y].vec;
                x1 = _mm256_min_epu16(x1, y1);
            }
        });
        *(__m256i *)(&randVec[x].items) = x1;
#else
        // 16 random 16-bit integers are realized by 4 64-bit random integers
        G->forNeighborsOf(x, [&](const node y) {
            if (distance[y] == distance[x] + 1)
                for (index i = 0; i < K; ++i)
                    randVec[K * x + i] = std::min(randVec[K * x + i], randVec[K * y + i]);
        });
#endif // AVX2_AVAILABLE

        if (distance[x] == 1) {
            sumOfMins[x] = 0;
            for (index j = 0; j < K; ++j) {
#ifdef AVX2_AVAILABLE
                sumOfMins[x] += randVec[x].items[j];
#else
                sumOfMins[x] += randVec[K * x + j];
#endif // AVX2_AVAILABLE
            }
            if (!sumOfMins[x])
                sumOfMins[x] = 1;
        }
    }

    G->forNodes([&](const node x) {
        if (distance[x] == 1) {
            float estimate =
                static_cast<float>(K) / (static_cast<float>(sumOfMins[x]) / maxInt16) - 1.f;
            if (estimate > bestEstimate) {
                v = x;
                bestEstimate = estimate;
            }
        }
    });

    assert(v != none);
    return v;
}

std::vector<node> GroupClosenessLocalSwaps::groupMaxCloseness() const {
    assureFinished();
    std::vector<node> maxGroup;
    maxGroup.reserve(group.size());

    for (const auto &entry : idxMap)
        maxGroup.push_back(entry.first);

    return maxGroup;
}

count GroupClosenessLocalSwaps::numberOfSwaps() const {
    assureFinished();
    return totalSwaps;
}

void GroupClosenessLocalSwaps::resetGamma(node x, index idx) {
    std::fill(gamma.begin() + group.size() * x, gamma.begin() + group.size() * (x + 1), false);
    gamma[group.size() * x + idx] = true;
}

bool GroupClosenessLocalSwaps::gammaIsU(node x, index idx) const {
    if (!gamma[group.size() * x + idx])
        return false;

    for (count i = 0; i < group.size(); ++i)
        if (i != idx && gamma[group.size() * x + i])
            return false;

    return true;
}

} // namespace NetworKit
