/*
 * PreferentialAttachmentIndex.cpp
 *
 *  Created on: 22.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "PreferentialAttachmentIndex.h"

namespace NetworKit {

PreferentialAttachmentIndex::PreferentialAttachmentIndex() {
}

PreferentialAttachmentIndex::PreferentialAttachmentIndex(const Graph& G) : LinkPredictor(G) {
}

double PreferentialAttachmentIndex::runImpl(node u, node v) {
  return G->degree(u) * G->degree(v);
}

} // namespace NetworKit
