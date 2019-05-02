#include <girgs/IntSort.h>
#include <girgs/ScopedTimer.h>
#include <girgs/Helper.h>


namespace girgs {


template<unsigned int D, typename EdgeCallback>
SpatialTree<D, EdgeCallback>::SpatialTree(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions, double alpha, EdgeCallback& edgeCallback, bool profile)
: m_EdgeCallback(edgeCallback)
, m_profile(profile)
, m_alpha(alpha)
, m_n(weights.size())
, m_w0(*std::min_element(weights.begin(), weights.end()))
, m_wn(*std::max_element(weights.begin(), weights.end()))
, m_W(std::accumulate(weights.begin(), weights.end(), 0.0))
, m_baseLevelConstant(static_cast<int>(std::log2(m_W/m_w0/m_w0))) // log2(W/w0^2)
, m_layers(static_cast<unsigned int>(floor(std::log2(m_wn/m_w0)))+1)
, m_levels(partitioningBaseLevel(0,0) + 1) // (log2(W/w0^2) - 2) / d
{
    assert(weights.size() == positions.size());
    assert(positions.size() > 0 && positions.front().size() == D);

    ScopedTimer timer("Preprocessing", profile);

    // determine which layer pairs to sample in which level
    m_layer_pairs.resize(m_levels);
    for (auto i = 0u; i < m_layers; ++i)
        for (auto j = 0u; j < m_layers; ++j)
            m_layer_pairs[partitioningBaseLevel(i, j)].emplace_back(i,j);

    // sort weights into exponentially growing layers
    {
        ScopedTimer timer("Build DS", profile);
        m_weight_layers = buildPartition(weights, positions);
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::generateEdges(int seed) {

    // one random generator and distribution for each thread
    const auto num_threads = omp_get_max_threads();
    m_gens.resize(num_threads);
    for (int thread = 0; thread < num_threads; thread++) {
        m_gens[thread].seed(seed >= 0 ? seed+thread : std::random_device()());
    } 

#ifndef NDEBUG
    // ensure that all node pairs are compared either type 1 or type 2
    m_type1_checks = 0;
    m_type2_checks = 0;
#endif // NDEBUG

    // sample all edges
	if (num_threads == 1) { 
        // sequential
		visitCellPair(0, 0, 0);
		assert(m_type1_checks + m_type2_checks == m_n*(m_n - 1ll));
		return;
    }

    // parallel see docs for visitCellPair_sequentialStart
    const auto first_parallel_level = static_cast<unsigned int>(std::ceil(std::log2(4.0*num_threads) / D));
    const auto parallel_cells = SpatialTreeCoordinateHelper<D>::numCellsInLevel(first_parallel_level);
    const auto first_parallel_cell = SpatialTreeCoordinateHelper<D>::firstCellOfLevel(first_parallel_level);

    // saw off recursion before "first_parallel_level" and save all calls that would be made
    auto parallel_calls = std::vector<std::vector<unsigned int>>(parallel_cells);
    visitCellPair_sequentialStart(0, 0, 0, first_parallel_level, parallel_calls);

    // do the collected calls in parallel
    #pragma omp parallel for schedule(static), num_threads(num_threads) // dynamic scheduling would be better but not reproducible
    for (int i = 0; i < parallel_cells; ++i) {
        auto current_cell = first_parallel_cell + i;
        for (auto each : parallel_calls[i])
            visitCellPair(current_cell, each, first_parallel_level);
    }

    assert(m_type1_checks + m_type2_checks == m_n*(m_n - 1ll));
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level) {
    if(!CoordinateHelper::touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
		if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    for(auto a = CoordinateHelper::firstChild(cellA); a<=CoordinateHelper::lastChild(cellA); ++a)
        for(auto b = cellA == cellB ? a : CoordinateHelper::firstChild(cellB); b<=CoordinateHelper::lastChild(cellB); ++b)
            visitCellPair(a, b, level+1);
}



template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
                                                   unsigned int first_parallel_level,
                                                   std::vector<std::vector<unsigned int>> &parallel_calls) {
    if(!CoordinateHelper::touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
		if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if (level == m_levels - 1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    for(auto a = CoordinateHelper::firstChild(cellA); a<=CoordinateHelper::lastChild(cellA); ++a)
        for(auto b = cellA == cellB ? a : CoordinateHelper::firstChild(cellB); b<=CoordinateHelper::lastChild(cellB); ++b){
            if(level+1 == first_parallel_level)
                parallel_calls[a-CoordinateHelper::firstCellOfLevel(first_parallel_level)].push_back(b);
            else
                visitCellPair_sequentialStart(a, b, level+1, first_parallel_level, parallel_calls);
        }
}



template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeI(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) == level || !CoordinateHelper::touching(cellA, cellB, level)); // in this case we were redirected from typeII with maxProb==1.0

    auto rangeA = m_weight_layers[i].cellIterators(cellA, level);
    auto rangeB = m_weight_layers[j].cellIterators(cellB, level);

    if (rangeA.first == rangeA.second || rangeB.first == rangeB.second)
        return;


#ifndef NDEBUG
    {
        const auto sizeV_i_A = std::distance(rangeA.first, rangeA.second);
        const auto sizeV_j_B = std::distance(rangeB.first, rangeB.second);

        #pragma omp atomic
        m_type1_checks += (cellA == cellB && i == j) ? sizeV_i_A * (sizeV_i_A - 1)  // all pairs in AxA without {v,v}
                                                     : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
    }
#endif // NDEBUG

    std::uniform_real_distribution<> dist;
    const auto threadId = omp_get_thread_num();

    const auto inThresholdMode = m_alpha == std::numeric_limits<double>::infinity();

    int kA = 0;
    for(auto pointerA = rangeA.first; pointerA != rangeA.second; ++kA, ++pointerA) {
        auto offset = (cellA == cellB && i==j) ? kA+1 : 0;
        for (auto pointerB = rangeB.first + offset; pointerB != rangeB.second; ++pointerB) {

            const auto& nodeInA = *pointerA;
            const auto& nodeInB = *pointerB;

			// pointer magic gives same results
			assert(nodeInA.index == m_weight_layers[i].kthPoint(cellA, level, kA).index);
			assert(nodeInB.index == m_weight_layers[j].kthPoint(cellB, level, std::distance(rangeB.first, pointerB)).index);

            // points are in correct cells
            assert(cellA - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(nodeInA.coord, level));
            assert(cellB - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(nodeInB.coord, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA.weight/m_w0)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB.weight/m_w0)));

            assert(nodeInA.index != nodeInB.index);
            const auto distance = nodeInA.distance(nodeInB);
            const auto w_term = nodeInA.weight*nodeInB.weight/m_W;
            const auto d_term = pow_to_the<D>(distance);

            if(inThresholdMode) {
                if(d_term < w_term)
                    m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
            } else {
                auto edge_prob = std::pow(w_term/d_term, m_alpha); // we don't need min with 1.0 here
                if(dist(m_gens[threadId]) < edge_prob)
                    m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
            }
        }
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeII(
        unsigned int cellA, unsigned int cellB, unsigned int level,
        unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) >= level);

    auto rangeA = m_weight_layers[i].cellIterators(cellA, level);
    auto rangeB = m_weight_layers[j].cellIterators(cellB, level);

    if (rangeA.first == rangeA.second || rangeB.first == rangeB.second)
        return;

    const auto sizeV_i_A = std::distance(rangeA.first, rangeA.second);
    const auto sizeV_j_B = std::distance(rangeB.first, rangeB.second);

    // get upper bound for probability
    const auto w_upper_bound = m_w0*(1<<(i+1)) * m_w0*(1<<(j+1)) / m_W;
    const auto cell_distance = CoordinateHelper::dist(cellA, cellB, level);
    const auto dist_lower_bound = pow_to_the<D>(cell_distance);
    const auto max_connection_prob = std::min(std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    const auto num_pairs = sizeV_i_A * sizeV_j_B;
    const auto expected_samples = num_pairs * max_connection_prob;

    // skipping over points is actually quite expensive as it messes up
    // branch predictions and prefetching. Hence low expected skip distances
    // it's cheapter to throw a coin each time!
    if (max_connection_prob > 0.2) {
        return sampleTypeI(cellA, cellB, level, i, j);
    }

#ifndef NDEBUG
    #pragma omp atomic
    m_type2_checks += 2 * sizeV_i_A * sizeV_j_B;
#endif // NDEBUG

    if(expected_samples < 1e-6)
        return;

    // init geometric distribution
    auto threadId = omp_get_thread_num();
    auto& gen = m_gens[threadId];
    auto geo = std::geometric_distribution<unsigned long long>(max_connection_prob);
    auto dist = std::uniform_real_distribution<>(0, max_connection_prob);

    for (auto r = geo(gen); r < num_pairs; r += 1 + geo(gen)) {
        // determine the r-th pair
        const Node<D>& nodeInA = rangeA.first[r%sizeV_i_A];
        const Node<D>& nodeInB = rangeB.first[r/sizeV_i_A];

        nodeInB.prefetch();
        nodeInA.prefetch();

        // points are in correct weight layer
        assert(i == static_cast<unsigned int>(std::log2(nodeInA.weight/m_w0)));
        assert(j == static_cast<unsigned int>(std::log2(nodeInB.weight/m_w0)));

        const auto rnd = dist(gen);

        // get actual connection probability
        const auto distance = nodeInA.distance(nodeInB);
        const auto w_term = nodeInA.weight*nodeInB.weight/m_W;
        const auto d_term = pow_to_the<D>(distance);
        const auto connection_prob = std::pow(w_term/d_term, m_alpha); // we don't need min with 1.0 here
        assert(w_term < w_upper_bound);
        assert(d_term >= dist_lower_bound);

        if(rnd < connection_prob) {
            m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
        }
    }
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::weightLayerTargetLevel(int layer) const {
    // -1 coz w0 is the upper bound for layer 0 in paper and our layers are shifted by -1
    auto result = std::max((m_baseLevelConstant - layer - 1) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct insertion level
        assert(0 <= layer && layer < m_layers);
        assert(0 <= result && result <= m_levels); // note the result may be one larger than the deepest level (hence the <= m_levels)
        auto volume_requested  = m_w0*m_w0*std::pow(2,layer+1)/m_W; // v(i) = w0*wi/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::partitioningBaseLevel(int layer1, int layer2) const {

    // we do the computation on signed ints but cast back after the max with 0
    // m_baseLevelConstant is just log(W/w0^2)
    auto result = std::max((m_baseLevelConstant - layer1 - layer2 - 2) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct comparison level
        assert(0 <= layer1 && layer1 < m_layers);
        assert(0 <= layer2 && layer2 < m_layers);
        auto volume_requested  = m_w0*std::pow(2,layer1+1) * m_w0*std::pow(2,layer2+1) / m_W; // v(i,j) = wi*wj/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}

template<unsigned int D, typename EdgeCallback>
std::vector<WeightLayer<D>> SpatialTree<D, EdgeCallback>::buildPartition(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions) {

    const auto n = weights.size();
    assert(positions.size() == n);

    auto weight_to_layer = [=] (double weight) {
        return std::log2(weight / m_w0);
    };

    const auto first_cell_of_layer = [&] {
        std::vector<unsigned int> first_cell_of_layer(m_layers + 1);
        unsigned int sum = 0;
        for (auto l = 0; l < m_layers; ++l) {
            first_cell_of_layer[l] = sum;
            sum += CoordinateHelper::numCellsInLevel(weightLayerTargetLevel(l));
        }
        first_cell_of_layer.back() = sum;
        return first_cell_of_layer;
    }();
    const auto max_cell_id = first_cell_of_layer.back();

    std::shared_ptr<Node<D>[]> points(new Node<D>[n]);
    // compute the cell a point belongs to
    {
        ScopedTimer timer("Classify points & precompute coordinates", m_profile);

        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            const auto layer = weight_to_layer(weights[i]);
            const auto level = weightLayerTargetLevel(layer);
            points[i] = Node<D>(positions[i], weights[i], i);
            points[i].cell_id = first_cell_of_layer[layer] + CoordinateHelper::cellForPoint(points[i].coord, level);
            assert(points[i].cell_id < max_cell_id);
        }
    }

    // Sort points by cell-ids
    {
        ScopedTimer timer("Sort points", m_profile);

        auto compare = [](const Node<D> &a, const Node<D> &b) { return a.cell_id < b.cell_id; };

        intsort::intsort(points, n, [](const Node<D> &p) { return p.cell_id; }, max_cell_id);
        //alternatively: std::sort(points.begin(), points.end(), compare);

        assert(std::is_sorted(points.get(), points.get() + n, compare));
    }


    // compute pointers into points
    constexpr auto gap_cell_indicator = std::numeric_limits<unsigned int>::max();
    std::shared_ptr<unsigned int[]> first_point_in_cell{new unsigned int[max_cell_id + 1]};
    std::fill_n(first_point_in_cell.get(), max_cell_id + 1, gap_cell_indicator);
    {
        ScopedTimer timer("Find first point in cell", m_profile);

        first_point_in_cell.get()[max_cell_id] = n;

        // First, we mark the begin of cells that actually contain points
        // and repair the gaps (i.e., empty cells) later. In the mean time,
        // the values of those gaps will remain at gap_cell_indicator.
        first_point_in_cell.get()[points.get()->cell_id] = 0;
        #pragma omp parallel for
        for (int i = 1; i < n; ++i) {
            if (points[i - 1].cell_id != points[i].cell_id) {
                first_point_in_cell.get()[points[i].cell_id] = i;
            }
        }

        // Now repair gaps: since first_point_in_cell shell contain
        // a prefix sum, we simply replace any "gap_cell_indicator"
        // with its nearest non-gap successor. In the main loop,
        // this is always the direct successors since we're iterating
        // from right to left.
        #pragma omp parallel
        {
            const auto threads = omp_get_num_threads();
            const auto rank = omp_get_thread_num();
            const auto chunk_size = (max_cell_id + threads - 1) / threads; // = ceil(max_cell_id / threads)

            // Fix right-most element (if gap) of each thread's elements by looking into chunk of next thread.
            // We do not need an end of array check, since it's guaranteed that the last element is n.
            // We're using on this very short code block to avoid UB even if we're only performing word-wise updates.
            #pragma omp single
            {
                for (int r = 0; r < threads; r++) {
                    const auto end = std::min(max_cell_id, chunk_size * (r + 1));
                    int first_non_invalid = end - 1;
                    while (first_point_in_cell.get()[first_non_invalid] == gap_cell_indicator)
                        first_non_invalid++;
                    first_point_in_cell.get()[end - 1] = first_point_in_cell.get()[first_non_invalid];
                }
            }

            const auto begin = std::min(max_cell_id, chunk_size * rank);

            auto i = std::min(max_cell_id, begin + chunk_size);
            while (i-- > begin) {
                first_point_in_cell.get()[i] = std::min(
                    first_point_in_cell.get()[i],
                    first_point_in_cell.get()[i + 1]);
            }
        }

        #ifndef NDEBUG
        {
            assert(points.get()[n-1].cell_id < max_cell_id);

            // assert that we have a prefix sum starting at 0 and ending in n
            assert(first_point_in_cell.get()[0] == 0);
            assert(first_point_in_cell.get()[max_cell_id] == n);
            assert(std::is_sorted(first_point_in_cell.get(), first_point_in_cell.get() + max_cell_id + 1));

            // check that each point is in its right cell (and that the cell boundaries are correct)
            for (auto cid = 0u; cid != max_cell_id; ++cid) {
                const auto begin = points.get() + first_point_in_cell.get()[cid];
                const auto end = points.get() + first_point_in_cell.get()[cid + 1];
                for (auto it = begin; it != end; ++it)
                    assert(it->cell_id == cid);
            }
        }
        #endif
    }

    // build spatial structure and find insertion level for each layer based on lower bound on radius for current and smallest layer
    std::vector<WeightLayer<D>> radius_layers;
    radius_layers.reserve(m_layers);
    {
        ScopedTimer timer("Build data structure", m_profile);
        for (auto layer = 0u; layer < m_layers; ++layer) {
            radius_layers.emplace_back(weightLayerTargetLevel(layer), points, first_point_in_cell, first_point_in_cell.get() + first_cell_of_layer[layer]);
        }
    }

    return radius_layers;
}


} // namespace girgs
