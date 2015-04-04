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
 *
 */
class PreferentialAttachmentIndex : public LinkPredictor {
private:
  /**
   *
   * @param u First node
   * @param v Second node
   * @return 
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* PREFERENTIALATTACHMENTINDEX_H_ */