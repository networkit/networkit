/*
 * LocalMover.h
 *
 *  Created on: 02.05.2013
 *      Author: cls
 */

#ifndef LOCALMOVER_H_
#define LOCALMOVER_H_

#include "Clusterer.h"

namespace NetworKit {

class LocalMover: public NetworKit::Clusterer {


public:

	class QualityObjective {
		virtual double getValue(node v) = 0;
	};

	class Modularity : public QualityObjective {
		virtual double getValue(node v);
	};

	class Coverage : public QualityObjective {
		virtual double getValue(node v);
	};


	LocalMover(QualityObjective& obj);

	virtual ~LocalMover();

	virtual Clustering run(Graph& G);

private:

	void move(node u, cluster C);

	QualityObjective* objective;
	Clustering* zeta;

};

} /* namespace NetworKit */
#endif /* LOCALMOVER_H_ */
