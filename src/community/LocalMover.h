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

	class Objective {

	public:

		/**
		 * Finds the best neighboring cluster for node  u.
		 */
		virtual cluster find(node u);

		virtual double eval(node u, cluster C);

		/**
		 * LocalMover notifies Objective to update data structures.
		 */
		virtual void onMove(node u, cluster C);


		class QualityObjective {
			virtual double getValue(node v) = 0;
		};

		class Modularity : public QualityObjective {
			virtual double getValue(node v);
		};

		class Coverage : public QualityObjective {
			virtual double getValue(node v);
		};

	};


//	class DeltaModularity : public Objective {
//
//	public:
//
//		virtual cluster find(node u);
//
//		virtual void onMove(node u, cluster C);
//
//	};

	class TerminationCriterion {
	public:

		virtual bool done();
	};

	LocalMover(Objective* obj, TerminationCriterion* crit);

	virtual ~LocalMover();

	virtual Clustering run(Graph& G);

private:

	void move(node u, cluster C);

	Objective* objective;
	TerminationCriterion* criterion;
	Clustering* zeta;

};

} /* namespace NetworKit */
#endif /* LOCALMOVER_H_ */
