#include <random>
#include <cmath>
#include <gtest/gtest.h>
#include <girgs/Helper.h>

template <typename T>
class Pow_test: public testing::Test {};

using test_types = ::testing::Types<
    std::integral_constant<std::size_t,0>,
    std::integral_constant<std::size_t,1>,
    std::integral_constant<std::size_t,2>,
    std::integral_constant<std::size_t,3>,
    std::integral_constant<std::size_t,5>,
    std::integral_constant<std::size_t,6>,
    std::integral_constant<std::size_t,7>,
    std::integral_constant<std::size_t,8>,
    std::integral_constant<std::size_t,9>
>;


TYPED_TEST_CASE(Pow_test, test_types);

TYPED_TEST(Pow_test, CrossTest) {
    static constexpr std::size_t D = TypeParam::value;

    std::mt19937_64 prng(D);
    std::uniform_real_distribution<double> dist(-1e2, 1e2);

    for(int i=0; i < 100000; ++i) {
        const auto x = dist(prng);
        const auto tested = girgs::pow_to_the<D>(x);

        auto ref = 1.0;
        for(int j=0; j < D; ++j) ref *= x;

        ASSERT_NEAR(tested, ref, fabs(ref * 1e-15));
    }
}
