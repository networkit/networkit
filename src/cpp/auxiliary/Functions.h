/*
 * Functions.h
 *
 *  Created on: 02.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <algorithm>

// INFO: see http://stackoverflow.com/questions/14199798/very-generic-argmax-function-in-c-wanted for discussion on argmax function

/**
 * Given
 * 		@param[in]	cont	a container
 * 		@param[in]	func	a function mapping an element from the container to a value
 * 	, return the element from the container for which the function maps to the maximum value.
 */
template<typename Container, typename F>
auto argmax(Container const & container, F && func) -> decltype(*std::begin(container))
{
    typedef decltype(*std::begin(container)) V;
    auto cmp = [&](V a, V b){ return func(a) < func(b); };
    return *std::max_element(std::begin(container), std::end(container), cmp);
}

/*
template<typename Iter, typename Function>
Iter argmax2(Iter begin, Iter end, Function f)
{
    typedef typename Iter::value_type T;
    return std::min_element(begin, end, [&value](const T& a, const T& b){
        return f(a) < f(b);
    });
}
*/



#endif /* FUNCTIONS_H_ */
