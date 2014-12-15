#ifndef FUNCTIONTRAITS_H
#define FUNCTIONTRAITS_H

#include <type_traits>
#include <tuple>

namespace Aux {

// Code by kennytm (auraHT Ltd.) 2011. See
// http://stackoverflow.com/questions/7943525/is-it-possible-to-figure-out-the-parameter-type-and-return-type-of-a-lambda
// Adapted by Michael Hamann <michael.hamann@kit.edu>
template <typename T>
struct FunctionTraits
	: public FunctionTraits<decltype(&T::operator())>
{};
// For generic types, directly use the result of the signature of its 'operator()'

template <typename ClassType, typename ReturnType, typename... Args>
struct FunctionTraits<ReturnType(ClassType::*)(Args...) const> {
	// we specialize for pointers to member function
	static constexpr std::size_t arity = sizeof...(Args);
	// arity is the number of arguments.

	using result_type = ReturnType;

	template <std::size_t i>
	struct arg;

	template <std::size_t i>
	struct arg {
			static_assert(i < arity, "error: invalid parameter index.");
			using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
			// the i-th argument is equivalent to the i-th tuple element of a tuple
			// composed of those arguments.
	};

};


} // namespace Aux

#endif // FUNCTIONTRAITS_H
