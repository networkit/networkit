#include <networkit/community/ParallelLeiden.hpp>

namespace NetworKit {
ParallelLeiden::ParallelLeiden(const Graph &graph, int iterations, bool randomize, double gamma)
    : CommunityDetectionAlgorithm(graph), gamma(gamma), numberOfIterations(iterations),
      random(randomize) {
    this->result = Partition(graph.numberOfNodes());
    this->result.allToSingletons();
}

void ParallelLeiden::run() {
    if (VECTOR_OVERSIZE < 1) {
        throw std::invalid_argument("VECTOR_OVERSIZE cant be smaller than 1");
    }
    auto totalTime = Aux::Timer();
    totalTime.start();
    do { // Leiden iteration
        INFO(numberOfIterations, " Leiden iteration(s) left");
        numberOfIterations--;
        changed = false;
        const Graph *currentGraph = G;
        Graph coarse;
        Partition refined;
        calculateVolumes(*currentGraph);
        do {
            handler.assureRunning();
            parallelMove(*currentGraph);
            // If each community consists of exactly one node we're done, i.e. when |V(G)| = |P|
            if (currentGraph->numberOfNodes() != result.numberOfSubsets()) {
                break;
            }
            handler.assureRunning();
            refined = parallelRefine(*currentGraph);
            handler.assureRunning();
            ParallelPartitionCoarsening ppc(*currentGraph, refined); // Aggregate graph
            ppc.run();
            auto temp = std::move(ppc.getCoarseGraph());
            auto map = std::move(ppc.getFineToCoarseNodeMapping());
            // Maintain Partition, add every coarse Node to the community its fine Nodes were in
            //  unlike in louvain, 2 coarse Nodes can belong to the same community
            Partition p(temp.numberOfNodes());
            p.setUpperBound(result.upperBound());
            currentGraph->parallelForNodes([&](node u) { p[map[u]] = result[u]; });

            mappings.emplace_back(std::move(map));
            result = std::move(p);
            std::swap(temp, coarse);
            currentGraph = &coarse;
        } while (true);
        flattenPartition();
        INFO("Leiden iteration done, took ", totalTime.elapsedTag(), "so far");
    } while (changed && numberOfIterations > 0);
    hasRun = true;
}

void ParallelLeiden::calculateVolumes(const Graph &graph) {
    auto timer = Aux::Timer();
    timer.start();
    // thread safe reduction. Avoid atomic calculation of total graph volume for unweighted graphs.
    // Vol(G) is then 2*|E|
    communityVolumes.clear();
    communityVolumes.resize(result.upperBound() + VECTOR_OVERSIZE);
    if (graph.isWeighted()) {
        std::vector<double> threadVolumes(omp_get_max_threads());
        graph.parallelForNodes([&](node a) {
            {
                edgeweight ew = graph.weightedDegree(a, true);
#pragma omp atomic
                communityVolumes[result[a]] += ew;
                threadVolumes[omp_get_thread_num()] += ew;
            }
        });
        for (const auto vol : threadVolumes) {
            inverseGraphVolume += vol;
        }
        inverseGraphVolume = 1 / inverseGraphVolume;
    } else {
        inverseGraphVolume = 1.0 / (2 * graph.numberOfEdges());
        graph.parallelForNodes([&](node a) {
            {
#pragma omp atomic
                communityVolumes[result[a]] += graph.weightedDegree(a, true);
            }
        });
    }
    TRACE("Calculating Volumes took " + timer.elapsedTag());
}

void ParallelLeiden::flattenPartition() {
    auto timer = Aux::Timer();
    timer.start();
    if (mappings.empty()) {
        return;
    }
    // Create a new partition with size |V(G)| (the fine/bigger Graph)
    Partition flattenedPartition(G->numberOfNodes());
    flattenedPartition.setUpperBound(result.upperBound());
    int i = mappings.size() - 1;
    std::vector<node> &lower = mappings[i--];
    while (i >= 0) {
        // iteratively "resolve" (i.e compose) mappings. Let "lower" be a mapping thats
        // below "higher" in the hierarchy (i.e. of a later aggregation)
        // If higher[index] = z and lower[z] = x then set higher[index] = x
        std::vector<node> &upper = mappings[i--];
        for (auto &idx : upper) {
            idx = lower[idx];
        }
        lower = upper;
    }
    G->parallelForNodes([&](node a) { flattenedPartition[a] = lower[a]; });
    flattenedPartition.compact(true);
    result = flattenedPartition;
    mappings.clear();
    TRACE("Flattening partition took " + timer.elapsedTag());
}

void ParallelLeiden::parallelMove(const Graph &graph) {
    DEBUG("Local Moving : ", graph.numberOfNodes(), " Nodes ");
    std::vector<count> moved(omp_get_max_threads(), 0);
    std::vector<count> totalNodesPerThread(omp_get_max_threads(), 0);
    std::atomic_int singleton(0);
    //^^^ Debug stuff

    // Only insert nodes to the queue when they're not already in it.
    std::vector<std::atomic_bool> inQueue(graph.upperNodeIdBound());
    std::queue<std::vector<node>> queue;
    std::mutex qlock;                      // queue lock
    std::condition_variable workAvailable; // waiting/notifying for new Nodes

    std::atomic_bool resize(false);
    std::atomic_int waitingForResize(0);
    std::atomic_int waitingForNodes(0);

    std::vector<int> order;
    int tshare;
    int tcount;
    uint64_t vectorSize = communityVolumes.capacity();
    std::atomic_int upperBound(result.upperBound());
#pragma omp parallel
    {
#pragma omp single
        {
            tcount = omp_get_num_threads();
            order.resize(tcount);
            for (int i = 0; i < tcount; i++) {
                order[i] = i;
            }
            if (random)
                std::shuffle(order.begin(), order.end(), Aux::Random::getURNG());
            tshare = 1 + graph.upperNodeIdBound() / tcount;
        }
        auto &mt = Aux::Random::getURNG();
        std::vector<node> currentNodes;
        currentNodes.reserve(tshare);
        std::vector<node> newNodes;
        newNodes.reserve(WORKING_SIZE);
        // cutWeight[Community] returns cut of Node to Community
        std::vector<double> cutWeights(communityVolumes.capacity());
        std::vector<index> pointers;
        int start = tshare * order[omp_get_thread_num()];
        int end = (1 + order[omp_get_thread_num()]) * tshare;

        for (int i = start; i < end; i++) {
            if (graph.hasNode(i)) {
                currentNodes.push_back(i);
                inQueue[i].store(true);
            }
        }
        if (random)
            std::shuffle(currentNodes.begin(), currentNodes.end(), mt);
#pragma omp barrier
        do {
            handler.assureRunning();
            for (node u : currentNodes) {
                // If a vector resize is needed, yield until done (This will probably never happen)
                // The causing thread will do so once waitingForResize == tcount (Amount of threads)
                if (resize) {
                    waitingForResize++;
                    while (resize) {
                        std::this_thread::yield();
                    }
                    waitingForResize--;
                }
                cutWeights.resize(vectorSize);
                assert(inQueue[u]);
                index currentCommunity = result[u];
                double maxDelta = std::numeric_limits<double>::lowest();
                index bestCommunity = none;
                double degree = 0;
                for (auto z : pointers) {
                    // Reset the clearlist : Set all cutweights to 0 and clear the pointer vector
                    cutWeights[z] = 0;
                }
                pointers.clear();

                graph.forNeighborsOf(u, [&](node neighbor, edgeweight ew) {
                    index neighborCommunity = result[neighbor];
                    if (cutWeights[neighborCommunity] == 0) {
                        pointers.push_back(neighborCommunity);
                    }
                    if (u == neighbor) {
                        degree += ew;
                    } else {
                        cutWeights[neighborCommunity] += ew;
                    }
                    degree += ew; // keep track of the nodes degree. Loops count twice
                });

                if (pointers.empty())
                    continue;

                // Determine Modularity delta for all neighbor communities
                for (auto community : pointers) {
                    // "Moving" a node to its current community is pointless
                    if (community != currentCommunity) {
                        double delta;
                        delta = modularityDelta(cutWeights[community], degree,
                                                communityVolumes[community]);
                        if (delta > maxDelta) {
                            maxDelta = delta;
                            bestCommunity = community;
                        }
                    }
                }
                double modThreshold = modularityThreshold(
                    cutWeights[currentCommunity], communityVolumes[currentCommunity], degree);

                if (0 > modThreshold || maxDelta > modThreshold) {
                    moved[omp_get_thread_num()]++;
                    if (0 > maxDelta) { // move node to empty community
                        singleton++;
                        bestCommunity = upperBound++;
                        if (bestCommunity >= communityVolumes.capacity()) {
                            // Wait until all other threads yielded, then increase vector size
                            // Chances are this will never happen. Ever. Seriously...
                            bool expected = false;
                            if (resize.compare_exchange_strong(expected, true)) {
                                vectorSize += VECTOR_OVERSIZE;
                                // this vector is thread-local so resizing is fine whenever.
                                cutWeights.resize(vectorSize);
                                while (waitingForResize < tcount - 1) {
                                    std::this_thread::yield();
                                }
                                // all other threads are yielding, so resize is fine
                                communityVolumes.resize(vectorSize);
                                expected = true;
                                resize.compare_exchange_strong(expected, false);
                            } else {
                                waitingForResize++;
                                while (resize) {
                                    std::this_thread::yield();
                                }
                                cutWeights.resize(vectorSize);
                                waitingForResize--;
                            }
                        }
                    }
                    result[u] = bestCommunity;
#pragma omp atomic
                    communityVolumes[bestCommunity] += degree;
#pragma omp atomic
                    communityVolumes[currentCommunity] -= degree;
                    changed = true;
                    bool expected = true;
                    inQueue[u].compare_exchange_strong(expected, false);
                    assert(expected);
                    graph.forNeighborsOf(u, [&](node neighbor) {
                        // Only add the node to the queue if it's not already
                        // in it, and it's not the Node we're currently moving
                        if (result[neighbor] != bestCommunity && neighbor != u) {
                            expected = false;
                            if (inQueue[neighbor].compare_exchange_strong(expected, true)) {
                                newNodes.push_back(neighbor);
                                // push new nodes to the queue in WORKING_SIZE steps
                                if (newNodes.size() == WORKING_SIZE) {
                                    qlock.lock();
                                    queue.emplace(std::move(newNodes));
                                    qlock.unlock();
                                    // Notify threads that new is available
                                    workAvailable.notify_all();
                                    newNodes.clear();
                                    newNodes.reserve(WORKING_SIZE);
                                }
                                assert(!expected);
                            }
                        }
                    });
                }
            }

            // queue check/wait
            // 3 cases : newnodes not empty -> continue, newnodes empty & queue not empty ->  pop
            // queue, both empty -> increment waiting, if waiting < #threads wait, else done notify
            // all
            totalNodesPerThread[omp_get_thread_num()] += currentNodes.size();
            if (!newNodes.empty()) {
                std::swap(currentNodes, newNodes);
                newNodes.clear();
                continue;
            }

            std::unique_lock<std::mutex> uniqueLock(qlock);
            if (!queue.empty()) {
                std::swap(currentNodes, queue.front());
                queue.pop();
            } else { // queue empty && newNodes empty
                waitingForNodes++;
                if (waitingForNodes < tcount) { // Not all nodes are done yet, wait for new work
                    waitingForResize++;
                    while (queue.empty() && waitingForNodes < tcount) {
                        workAvailable.wait(uniqueLock);
                    }
                    if (waitingForNodes < tcount) {
                        // Notified and not all done means there's new work
                        std::swap(currentNodes, queue.front());
                        queue.pop();
                        waitingForNodes--;
                        waitingForResize--;
                        continue;
                    }
                }
                uniqueLock.unlock(); // Notified and all done, stop.
                workAvailable.notify_all();
                break;
            }
        } while (true);
        TRACE("Thread ", omp_get_thread_num(), " worked ",
              totalNodesPerThread[omp_get_thread_num()], "Nodes and moved ",
              moved[omp_get_thread_num()]);
    }
    result.setUpperBound(upperBound);
    assert(queue.empty());
    assert(waitingForNodes == tcount);
    if (Aux::Log::isLogLevelEnabled(Aux::Log::LogLevel::DEBUG)) {
        count totalMoved = std::accumulate(moved.begin(), moved.end(), (count)0);
        count totalWorked =
            std::accumulate(totalNodesPerThread.begin(), totalNodesPerThread.end(), (count)0);
        tlx::unused(totalMoved); // get around unused variable
        tlx::unused(totalWorked);
        DEBUG("Total worked: ", totalWorked, " Total moved: ", totalMoved,
              " moved to singleton community: ", singleton);
    }
}

Partition ParallelLeiden::parallelRefine(const Graph &graph) {
    Partition refined(graph.numberOfNodes());
    refined.allToSingletons();
    DEBUG("Starting refinement with ", result.numberOfSubsets(), " partitions");
    std::vector<uint_fast8_t> singleton(refined.upperBound(), true);
    std::vector<double> cutCtoSminusC(refined.upperBound());
    std::vector<double> refinedVolumes(refined.upperBound()); // Community Volumes P_refined
    std::vector<std::mutex> locks(refined.upperBound());
    std::vector<node> nodes(graph.upperNodeIdBound(), none);
#pragma omp parallel
    {
        std::vector<index> neighComms;
        // Keeps track of relevant Neighbor communities. Needed to reset the clearlist fast
        std::vector<double> cutWeights(refined.upperBound()); // cut from Node to Communities
        auto &mt = Aux::Random::getURNG();
#pragma omp for
        for (omp_index u = 0; u < static_cast<omp_index>(graph.upperNodeIdBound()); u++) {
            if (graph.hasNode(u)) {
                nodes[u] = u;
                graph.forNeighborsOf(u, [&](node neighbor, edgeweight ew) {
                    if (u != neighbor) {
                        if (result[neighbor] == result[u]) {
                            // Cut to communities in the refined partition that
                            // are in the same community in the original partition
                            cutCtoSminusC[u] += ew;
                        }
                    } else {
                        refinedVolumes[u] += ew;
                    }
                    refinedVolumes[u] += ew;
                });
            }
        }
        if (random) {
            int share = graph.upperNodeIdBound() / omp_get_num_threads();
            int start = omp_get_thread_num() * share;
            int end = (omp_get_thread_num() + 1) * share - 1;
            if (omp_get_thread_num() == omp_get_num_threads() - 1)
                end = nodes.size() - 1;
            if (start != end && end > start)
                std::shuffle(nodes.begin() + start, nodes.begin() + end, mt);
#pragma omp barrier
        }
        handler.assureRunning();
#pragma omp for schedule(dynamic, WORKING_SIZE)
        for (omp_index i = 0; i < static_cast<omp_index>(nodes.size()); i++) {
            node u = nodes[i];
            if (u == none || !singleton[u]) { // only consider singletons
                continue;
            }
            index S = result[u];                // Node's community ID in the original partition (S)
            for (auto neighComm : neighComms) { // Reset the clearlist : Set all cutweights to 0
                if (neighComm != none)
                    cutWeights[neighComm] = 0;
            }

            neighComms.clear();

            std::vector<node> criticalNodes;
            // Nodes whose community ID equals their Node ID. These may be singletons that can
            // affect the cut which we need to update later
            double degree = 0;

            graph.forNeighborsOf(u, [&](node neighbor, edgeweight ew) { // Calculate degree and cut
                degree += ew;
                if (neighbor != u) {
                    if (S == result[neighbor]) {
                        index z = refined[neighbor];
                        if (z == neighbor) {
                            criticalNodes.push_back(neighbor);
                        }
                        // We don't need to remember the weight of that edge since it's
                        // already saved in cutWeights
                        if (cutWeights[z] == 0)
                            neighComms.push_back(z); // Keep track of neighbor communities
                        cutWeights[z] += ew;
                    }
                } else {
                    degree += ew;
                }
            });
            if (cutCtoSminusC[u] < this->gamma * degree * (communityVolumes[S] - degree)
                                       * inverseGraphVolume) { // R-Set Condition
                continue;
            }

            if (cutWeights[u] != 0) { // Node has been moved -> not a singleton anymore. Stop.
                continue;
            }

            double delta;
            index bestC = none;
            double bestDelta = std::numeric_limits<double>::lowest();
            int idx;
            // Determine Community that yields highest modularity delta
            auto bestCommunity = [&] {
                for (unsigned int j = 0; j < neighComms.size(); j++) {
                    index C = neighComms[j];
                    if (C == none) {
                        continue;
                    }
                    delta = modularityDelta(cutWeights[C], degree, refinedVolumes[C]);

                    if (delta < 0) { // modThreshold is 0, since cutw(v,C-) = 0 and volw(C-) = 0
                        continue;
                    }

                    auto volC = refinedVolumes[C];
                    if (delta > bestDelta
                        && cutCtoSminusC[C] >= this->gamma * volC * (communityVolumes[S] - volC)
                                                   * inverseGraphVolume) { // T-Set Condition
                        bestDelta = delta;
                        bestC = C;
                        idx = j;
                    }
                }
            };
            // To update cut values in case neighbors of this node have moved in the meantime
            auto updateCut = [&] {
                for (node &neighbor : criticalNodes) {
                    if (neighbor != none) {
                        index neighborCommunity = refined[neighbor];
                        if (neighborCommunity != neighbor) {
                            if (cutWeights[neighborCommunity] == 0) {
                                // remember to clear the vector, this community was not saved
                                // initially since the neighbor moved to it later
                                neighComms.push_back(neighborCommunity);
                            }
                            // cutWeights[Neighbor] is the weight of the edge between Node and
                            // Neighbor, since Neighbor was a singleton
                            cutWeights[neighborCommunity] += cutWeights[neighbor];
                            // Clear cutWeights entry beforehand, so we can "erase" bestC
                            // from the pointers vector by replacing it with "none"
                            cutWeights[neighbor] = 0;
                            neighbor = none;
                        }
                    }
                }
            };
            bestCommunity();
            if (bestC == none) {
                continue;
            }
            lockLowerFirst(u, bestC, locks); // avoid deadlocks
            // If this node is no longer a singleton, stop.
            if (singleton[u]) {
                // Target community still contains its "host" node? If not,this community is
                // now empty, choose a new one.
                while (bestC != none && refined[bestC] != bestC) {
                    locks[bestC].unlock();
                    // This makes sure it won't be considered in the next bestCommunity() call
                    neighComms[idx] = none;
                    bestC = none;
                    bestDelta = std::numeric_limits<double>::lowest();
                    updateCut();
                    bestCommunity();
                    if (bestC != none) {
                        if (!locks[bestC].try_lock()) {
                            if (u < bestC) {
                                locks[bestC].lock();
                            } else {
                                // temporarily release lock on current Node to avoid deadlocks
                                locks[u].unlock();
                                lockLowerFirst(u, bestC, locks);
                            }
                            if (!singleton[u]) {
                                locks[u].unlock();
                                locks[bestC].unlock();
                                continue;
                            }
                        }
                    }
                }
                if (bestC == none) {
                    locks[u].unlock(); // bestC was already unlocked in the loop above
                    continue;
                }
                singleton[bestC] = false;
                refined[u] = bestC;
                refinedVolumes[bestC] += degree;
                updateCut();
                cutCtoSminusC[bestC] += cutCtoSminusC[u] - 2 * cutWeights[bestC];
            }
            locks[bestC].unlock();
            locks[u].unlock();
        }
    }

    DEBUG("Ending refinement with ", refined.numberOfSubsets(), " partitions");
    return refined;
}
} // namespace NetworKit
