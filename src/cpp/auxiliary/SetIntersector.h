/*
 * SetIntersector.h
 *
 *  Created on: 14.05.2014
 *      Author: Henning
 */

#ifndef SETINTERSECTOR_H_
#define SETINTERSECTOR_H_

#include <vector>
#include <set>

namespace Aux {

template<class T>
class SetIntersector {
public:
	SetIntersector(T upperBound);

	/**
	 * @return Intersection of sets provided in @a A and @a B.
	 */
	std::set<T> intersect(const std::vector<T>& A, const std::vector<T>& B);

private:
	std::vector<bool> bv;
	uint64_t n;
};

} // namespace Aux


template<class T>
inline Aux::SetIntersector<T>::SetIntersector(T upperBound): n(upperBound) {
	bv.resize(n, false);

}

template<class T>
inline std::set<T> Aux::SetIntersector<T>::intersect(const std::vector<T>& A,
		const std::vector<T>& B)
{
	const std::vector<T>& smaller = (A.size() <= B.size()) ? A : B;
	const std::vector<T>& larger = (A.size() <= B.size()) ? B : A;

	for (auto entry: smaller) {
		bv[entry] = true;
	}

	std::set<T> result;
	for (auto entry: larger) {
		if (bv[entry]) {
			result.insert(entry);
		}
	}

	for (auto entry: smaller) {
		bv[entry] = false;
	}

	return result;
}

#endif /* SETINTERSECTOR_H_ */
