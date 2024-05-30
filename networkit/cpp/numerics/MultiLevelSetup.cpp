/*
 * MultiLevelSetup.cpp
 *
 *  Created on: 23.04.2024
 *      Author: Lukas
 */

#include <networkit/numerics/LAMG/MultiLevelSetup.hpp>

namespace NetworKit {

template class MultiLevelSetup<CSRMatrix>;
template class MultiLevelSetup<DenseMatrix>;
template class MultiLevelSetup<DynamicMatrix>;

} // namespace NetworKit
