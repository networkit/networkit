/*
 * PreferentialAttachmentIndex.h
 *
 *  Created on: 22.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef PREFERENTIALATTACHMENTINDEX_H_
#define PREFERENTIALATTACHMENTINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Preferential Attachment Index. 
 * The run-method simply calculates the product of the number of nodes in the neighborhoods
 * regarding the given nodes.
 */
class PreferentialAttachmentIndex : public LinkPredictor {
private:
  /**
   * Returns the product of the cardinalities of the neighborhoods regarding @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the product of the cardinalities of the neighborhoods regarding @a u and @a v
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* PREFERENTIALATTACHMENTINDEX_H_ */