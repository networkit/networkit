/*
 * Lamg.cpp
 *
 *  Created on: 02.05.2024
 */

#include <networkit/numerics/LAMG/Lamg.hpp>

namespace NetworKit {

template class Lamg<CSRMatrix>;
template class Lamg<DenseMatrix>;
template class Lamg<DynamicMatrix>;

} // namespace NetworKit
