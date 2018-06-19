/*
 * CurveballGlobalTradeSequence.h
 *
 *  Created on: 23.05.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef CURVEBALL_GLOBAL_TRADE_SEQUENCE_H_
#define CURVEBALL_GLOBAL_TRADE_SEQUENCE_H_

#include <cmath>
#include <tuple>
#include <type_traits>
#include <random>
#include <cassert>
#include <iostream>

#include "../Globals.h"
#include "../auxiliary/Log.h"

namespace NetworKit {
namespace CurveballDetails {

/**
 * Computes a bijection f:[p]->[p], x -> (a*x+b) mod p where p is chosen
 * as the smallest p >= n and p prime. a and b are drawn unif at random.
 * @tparam T
 */
template <typename T>
class LinearCongruentialMap {
    static_assert(!std::is_signed<T>::value,
                  "LinearCongruentialMap requires unsigned types");

    using signed_value_type = typename std::make_signed<T>::type;
    using signed_tuple = std::tuple<signed_value_type, signed_value_type, signed_value_type>;

public:
    static constexpr bool has_invert = true;
    using value_type = T;

    LinearCongruentialMap() {}

    explicit LinearCongruentialMap(value_type n, unsigned seed = 0) :
        n_(n), p_(compute_next_prime_(n))
    {
        std::mt19937_64 prng(seed * n + n);
        sample_parameters(prng);
    }

    LinearCongruentialMap(value_type n, std::mt19937_64 & prng) :
        n_(n), p_(compute_next_prime_(n))
    {
        sample_parameters(prng);
    }


    LinearCongruentialMap(value_type n, value_type a, value_type b) :
        n_(n), p_(compute_next_prime_(n)), a_(a),
        ainv_(static_cast<value_type>((std::get<1>(gcd_extended_(a_, p_)) + p_) % p_)),
        b_(b)
    {}

    LinearCongruentialMap(const LinearCongruentialMap&) = default;
    LinearCongruentialMap(LinearCongruentialMap&&) = default;
    LinearCongruentialMap& operator=(const LinearCongruentialMap&) = default;
    LinearCongruentialMap& operator=(LinearCongruentialMap&&) = default;


    //! Hashes [n] to [p] with (a*x+b) mod p
    value_type hash(value_type n) const {
        return (a_ * n + b_) % p_;
    }

    //! Alias to hash(n)
    value_type operator() (value_type n) const {
        return hash(n);
    }


    //! invert(hash(x)) == x
    value_type invert(value_type y) const {
        return ainv_ * (y + p_ - b_) % p_;
    }

    //! Is true if there exists no x in [n], s.t. h(x) == y,
    //! i.e. it can be used to determine if the map maps to y
    bool is_gap(value_type y) const {
        return invert(y) >= n_;
    }


    //! randomly samples parameters a and b
    void sample_parameters(std::mt19937_64& prng) {
        {
            const value_type max_a = std::numeric_limits<value_type>::max() / n_ - 1;
            if (max_a < p_) {
                std::cerr << "WARNING: Reduce randomness of hash function to avoid integer precision issues\n";
            }


            std::uniform_int_distribution<value_type> distr(1, std::min<value_type>(p_,  max_a) - 1);
            a_ = distr(prng);
            ainv_ = static_cast<value_type>( (std::get<1>(gcd_extended_(a_, p_)) + p_) % p_ );
        }
        {
            std::uniform_int_distribution<value_type> distr(0, p_ - 1);
            b_ = distr(prng);
        }
    }


    //! Sets parameters a = 1 and b = 0
    void set_as_identity() {
        a_ = 1;
        b_ = 0;
        ainv_ = 1;
    }

    value_type param_a() const {return a_;}
    value_type param_ainv() const {return ainv_;}
    value_type param_b() const {return b_;}
    value_type param_p() const {return p_;}

private:
    value_type n_; //< number of elements to be mapped
    value_type p_; //< size of the universe chosen to be >= n

    value_type a_; //< multiplicative parameter in hash
    value_type ainv_; //< multiplicative invert

    value_type b_; //< additive parameter in hash

    value_type compute_next_prime_(value_type n) const {
        auto is_prime = [](value_type n) {
            if (n <= 3) return true;
            if (0 == n % 2 || 0 == n % 3) return false;

            const value_type sqrt = static_cast<value_type>(std::sqrt(n) + 2);
            for (value_type i = 5; i < sqrt; i += 6)
                if (0 == n % i) return false;

            return true;
        };

        while (!is_prime(n))
            n++;

        return n;
    }

    // extended euclidian algorithm with 1 = gcd(a, p) = a*s + t*p mod p = a*s
    // --> s = 1/a
    signed_tuple gcd_extended_(const signed_value_type a, const signed_value_type b) const
    {
        if (a == 0)
            return signed_tuple{b, 0, 1};

        const value_type div = b / a;
        const value_type rem = b % a;

        auto tmp = gcd_extended_(rem, a);
        value_type x = std::get<2>(tmp) - div * std::get<1>(tmp);
        auto result = std::make_tuple(std::get<0>(tmp), x, std::get<1>(tmp));

        return result;
    }
};


/**
 * Computes a bijection f:[p]->[p], x -> (a*x+b) mod p where p is chosen
 * as the smallest p >= n and p prime. a and b are drawn unif at random.
 * @tparam T
 */
template <typename T>
class FixedLinearCongruentialMap {
    static_assert(!std::is_signed<T>::value,
                  "LinearCongruentialMap requires unsigned types");

    using signed_value_type = typename std::make_signed<T>::type;
    using signed_tuple = std::tuple<signed_value_type, signed_value_type, signed_value_type>;

public:
    static constexpr bool has_invert = true;
    using value_type = T;

    FixedLinearCongruentialMap() {}

    explicit FixedLinearCongruentialMap(value_type n, unsigned seed = 0) :
        n_(n)
    {
        std::mt19937_64 prng(seed * n + n);
        sample_parameters(prng);
    }

    FixedLinearCongruentialMap(value_type n, std::mt19937_64 & prng) :
        n_(n)
    {
        sample_parameters(prng);
    }


    FixedLinearCongruentialMap(value_type n, value_type a, value_type b) :
        n_(n), a_(a),
        ainv_(static_cast<value_type>((std::get<1>(gcd_extended_(a_, p_)) + p_) % p_)),
        b_(b)
    {}

    FixedLinearCongruentialMap(const FixedLinearCongruentialMap&) = default;
    FixedLinearCongruentialMap(FixedLinearCongruentialMap&&) = default;
    FixedLinearCongruentialMap& operator=(const FixedLinearCongruentialMap&) = default;
    FixedLinearCongruentialMap& operator=(FixedLinearCongruentialMap&&) = default;


    //! Hashes [n] to [p] with (a*x+b) mod p
    value_type hash(value_type n) const {
        return (a_ * n + b_) % p_;
    }

    //! Alias to hash(n)
    value_type operator() (value_type n) const {
        return hash(n);
    }


    //! invert(hash(x)) == x
    value_type invert(value_type y) const {
        return ainv_ * (y + p_ - b_) % p_;
    }

    //! Is true if there exists no x in [n], s.t. h(x) == y,
    //! i.e. it can be used to determine if the map maps to y
    bool is_gap(value_type y) const {
        return invert(y) >= n_;
    }


    //! randomly samples parameters a and b
    void sample_parameters(std::mt19937_64& prrg) {
        {
            std::uniform_int_distribution<value_type> distr(1, p_ - 1);
            a_ = distr(prrg);
            ainv_ = static_cast<value_type>( (std::get<1>(gcd_extended_(a_, p_)) + p_) % p_ );
        }
        {
            std::uniform_int_distribution<value_type> distr(0, p_ - 1);
            b_ = distr(prrg);
        }
    }


    //! Sets parameters a = 1 and b = 0
    void set_as_identity() {
        a_ = 1;
        b_ = 0;
        ainv_ = 1;
    }

    value_type param_a() const {return a_;}
    value_type param_ainv() const {return ainv_;}
    value_type param_b() const {return b_;}
    value_type param_p() const {return p_;}

private:
    value_type n_; //< number of elements to be mapped
    static constexpr value_type p_ = (1 || sizeof(value_type) == 4) ? 2147483647 : 2305843009213693951;

    value_type a_; //< multiplicative parameter in hash
    value_type ainv_; //< multiplicative invert

    value_type b_; //< additive parameter in hash

    // extended euclidian algorithm with 1 = gcd(a, p) = a*s + t*p mod p = a*s
    // --> s = 1/a
    signed_tuple gcd_extended_(const signed_value_type a, const signed_value_type b) const
    {
        if (a == 0)
            return signed_tuple{b, 0, 1};

        const value_type div = b / a;
        const value_type rem = b % a;

        auto tmp = gcd_extended_(rem, a);
        value_type x = std::get<2>(tmp) - div * std::get<1>(tmp);
        auto result = std::make_tuple(std::get<0>(tmp), x, std::get<1>(tmp));

        return result;
    }
};

/**
 * Handles a sequence of num_global_trades hash functions each with its own
 * set of parameters. It prominently exposes the "current" hash function
 * and its successors. Switching the current hash function is possible with
 * using the switch_to_round method.
 *
 * hash, invert  access the "current" hash function
 * hash_next     accesses the successor of "current" which is the identity if
 *               there are no further hash functions
 *
 * @tparam Default constructable, copy-able, requires methods
 *         sample_parameters, set_as_identity, hash and (optinally) invert
 */
template <typename Hash>
class GlobalTradeSequence {
public:
    using value_type = typename Hash::value_type;

    GlobalTradeSequence(node num_nodes, size_t num_global_trades, std::mt19937_64& prng) {
        hash_functors_.reserve(num_global_trades);
        hash_functors_.push_back(Hash{num_nodes, prng});

        for(size_t i=0; i < num_global_trades; i++) {
            // we copy the hash function, rather than constructing a new one,
            // to avoid repeated computations, such as the next larger prime
            // number
            hash_functors_.push_back(hash_functors_.back());
            hash_functors_.back().sample_parameters(prng);
        }

        switch_to_round(0);
    }

    void switch_to_round(size_t round) {
        assert(round < hash_functors_.size());
        current_ = hash_functors_[round];
        if (round + 1 == hash_functors_.size()) {
            next_ = current_;
            next_.set_as_identity();
        } else {
            next_ = hash_functors_[round + 1];
        }
    }

    value_type hash(node u) const {
        num_hashed_++;
        return current_.hash(u);
    }

    node invert(value_type u) const {
        num_inverted_++;
        return current_.invert(u);
    }

    value_type hash_next(node u) const {
        num_hashed_++;
        return next_.hash(u);
    }

    size_t number_of_rounds() const {
        return hash_functors_.size();
    }

    ~GlobalTradeSequence() {
        INFO("Hashed: ", num_hashed_, " Inverted: ", num_inverted_);
    }

private:
    std::vector<Hash> hash_functors_; // Todo: make tlx::simple_vector
    Hash current_;
    Hash next_;

    mutable count num_hashed_{0};
    mutable count num_inverted_{0};
};


} // curveball_details
} // namespace NetworKit


#endif // !CURVEBALL_GLOBAL_TRADE_SEQUENCE_H_
