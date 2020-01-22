/*
 * PreferentialAttachmentIndex.cpp
 *
 *  Created on: 22.03.2015
 *      Author: Kolja Esders
 */

#include <networkit/linkprediction/PreferentialAttachmentIndex.hpp>

namespace NetworKit {

double PreferentialAttachmentIndex::runImpl(node u, node v) {
  return G->degree(u) * G->degree(v);
}

} // namespace NetworKit
