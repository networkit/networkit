/*
 * LevelHierarchy.cpp
 *
 *  Created on: 23.04.2024
 *      Author: Lukas
 */

#include <networkit/numerics/LAMG/LevelHierarchy.hpp>

namespace NetworKit {

template class LevelHierarchy<CSRMatrix>;
template class LevelHierarchy<DenseMatrix>;
template class LevelHierarchy<DynamicMatrix>;

} // namespace NetworKit
