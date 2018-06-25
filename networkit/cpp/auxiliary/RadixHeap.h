// This file contains several Header files from https://github.com/tlx/tlx

 /*******************************************************************************
 * tlx/container/radixheap.hpp
 *
 * Part of tlx - http://panthema.net/tlx
 *
 * Copyright (C) 2018 Manuel Penschuck <tlx@manuel.jetzt>
 *
 * All rights reserved. Published under the Boost Software License, Version 1.0
 ******************************************************************************/

#ifndef TLX_CONTAINER_RADIXHEAP_HEADER
#define TLX_CONTAINER_RADIXHEAP_HEADER

#include <array>
#include <cassert>
#include <limits>
#include <type_traits>
#include <utility>

#include <vector>
#include <deque>

#include <cstdint>

#ifdef _MSC_VER
#include <intrin.h>
#endif


namespace Aux {
namespace radixheap_detail {
/*******************************************************************************
 * tlx/meta/log2.hpp
 *
 * Template Metaprogramming Tools (from the Generative Programming book Krysztof
 * Czarnecki, Ulrich Eisenecker)
 *
 * Part of tlx - http://panthema.net/tlx
 *
 * Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 * Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 * Copyright (C) 2010-2017 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the Boost Software License, Version 1.0
 ******************************************************************************/

#ifndef TLX_META_LOG2_HEADER
#define TLX_META_LOG2_HEADER


namespace tlx {

//! \addtogroup tlx_meta
//! \{

/******************************************************************************/
// Log2Floor<Value>::value

template <uint64_t Input>
class Log2Floor
{
public:
    enum {
        value = Log2Floor<Input / 2>::value + 1
    };
};

template <>
class Log2Floor<1>
{
public:
    enum { value = 0 };
};

template <>
class Log2Floor<0>
{
public:
    enum { value = 0 };
};

/******************************************************************************/
// Log2<Value>::floor and Log2<Value>::ceil

template <uint64_t Input>
class Log2
{
public:
    enum {
        floor = Log2Floor<Input>::value,
        ceil = Log2Floor<Input - 1>::value + 1
    };
};

template <>
class Log2<1>
{
public:
    enum {
        floor = 0,
        ceil = 0
    };
};

template <>
class Log2<0>
{
public:
    enum {
        floor = 0,
        ceil = 0
    };
};

//! \}

} // namespace tlx

#endif // !TLX_META_LOG2_HEADER

/******************************************************************************/


/*******************************************************************************
 * tlx/math/clz.hpp
 *
 * clz() count leading zeros - mainly for portability.
 *
 * Part of tlx - http://panthema.net/tlx
 *
 * Copyright (C) 2017 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the Boost Software License, Version 1.0
 ******************************************************************************/

#ifndef TLX_MATH_CLZ_HEADER
#define TLX_MATH_CLZ_HEADER


namespace tlx {

//! \addtogroup tlx_math
//! \{

/******************************************************************************/
// clz() - count leading zeros

//! clz (count leading zeros) - generic implementation
template <typename Integral>
static inline unsigned clz_template(Integral x) {
    if (x == 0) return 8 * sizeof(x);
    unsigned r = 0;
    while ((x & (static_cast<Integral>(1) << (8 * sizeof(x) - 1))) == 0)
        x <<= 1, ++r;
    return r;
}

/******************************************************************************/

template <typename Integral>
inline unsigned clz(Integral x);

#if defined(__GNUC__) || defined(__clang__)

//! clz (count leading zeros)
template <>
inline unsigned clz<unsigned>(unsigned i) {
    return static_cast<unsigned>(__builtin_clz(i));
}

//! clz (count leading zeros)
template <>
inline unsigned clz<int>(int i) {
    return clz(static_cast<unsigned>(i));
}

//! clz (count leading zeros)
template <>
inline unsigned clz<unsigned long>(unsigned long i) {
    return static_cast<unsigned>(__builtin_clzl(i));
}

//! clz (count leading zeros)
template <>
inline unsigned clz<long>(long i) {
    return clz(static_cast<unsigned long>(i));
}

//! clz (count leading zeros)
template <>
inline unsigned clz<unsigned long long>(unsigned long long i) {
    return static_cast<unsigned>(__builtin_clzll(i));
}

//! clz (count leading zeros)
template <>
inline unsigned clz<long long>(long long i) {
    return clz(static_cast<unsigned long long>(i));
}

#elif defined(_MSC_VER)

//! clz (count leading zeros)
template <typename Integral>
inline unsigned clz<unsigned>(Integral i) {
    unsigned long leading_zeros = 0;
    if (sizeof(i) > 4) {
        if (_BitScanReverse64(&leading_zeros, i))
            return 63 - leading_zeros;
        else
            return 8 * sizeof(i);
    }
    else {
        if (_BitScanReverse(&leading_zeros, static_cast<unsigned>(i)))
            return 31 - leading_zeros;
        else
            return 8 * sizeof(i);
    }
}

#else

//! clz (count leading zeros)
template <>
inline unsigned clz<int>(int i) {
    return clz_template(i);
}

//! clz (count leading zeros)
template <>
inline unsigned clz<unsigned>(unsigned i) {
    return clz_template(i);
}

//! clz (count leading zeros)
template <>
inline unsigned clz<long>(long i) {
    return clz_template(i);
}

//! clz (count leading zeros)
template <>
inline unsigned clz<unsigned long>(unsigned long i) {
    return clz_template(i);
}

//! clz (count leading zeros)
template <>
inline unsigned clz<long long>(long long i) {
    return clz_template(i);
}

//! clz (count leading zeros)
template <>
inline unsigned clz<unsigned long long>(unsigned long long i) {
    return clz_template(i);
}

#endif

//! \}

} // namespace tlx

#endif // !TLX_MATH_CLZ_HEADER

/******************************************************************************/

/*******************************************************************************
 * tlx/math/ffs.hpp
 *
 * ffs() find first set bit in integer - mainly for portability as ffs() is a
 * glibc extension and not available on Visual Studio.
 *
 * Part of tlx - http://panthema.net/tlx
 *
 * Copyright (C) 2007-2017 Timo Bingmann <tb@panthema.net>
 *
 * All rights reserved. Published under the Boost Software License, Version 1.0
 ******************************************************************************/

#ifndef TLX_MATH_FFS_HEADER
#define TLX_MATH_FFS_HEADER

namespace tlx {

//! \addtogroup tlx_math
//! \{

/******************************************************************************/
// ffs() - find first set bit in integer

//! ffs (find first set bit) - generic implementation
template <typename Integral>
static inline unsigned ffs_template(Integral x) {
    if (x == 0) return 0u;
    unsigned r = 1;
    while ((x & 1) == 0)
        x >>= 1, ++r;
    return r;
}

/******************************************************************************/

#if defined(__GNUC__) || defined(__clang__)

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(int i) {
    return static_cast<unsigned>(__builtin_ffs(i));
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(unsigned i) {
    return ffs(static_cast<int>(i));
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(long i) {
    return static_cast<unsigned>(__builtin_ffsl(i));
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(unsigned long i) {
    return ffs(static_cast<long>(i));
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(long long i) {
    return static_cast<unsigned>(__builtin_ffsll(i));
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(unsigned long long i) {
    return ffs(static_cast<long long>(i));
}

#else

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(int i) {
    return ffs_template(i);
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(unsigned int i) {
    return ffs_template(i);
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(long i) {
    return ffs_template(i);
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(unsigned long i) {
    return ffs_template(i);
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(long long i) {
    return ffs_template(i);
}

//! find first set bit in integer, or zero if none are set.
static inline
unsigned ffs(unsigned long long i) {
    return ffs_template(i);
}

#endif

//! \}

} // namespace tlx

#endif // !TLX_MATH_FFS_HEADER

/******************************************************************************/



/*!
 * Compute the rank of an integer x (i.e. the number of elements smaller than x
 * that are representable using type Int) and vice versa.
 * If Int is an unsigned integral type, all computations yield identity.
 * If Int is a signed integrals, the smallest (negative) number is mapped to
 * rank zero, the next larger value to one and so on.
 *
 * The implementation assumes negative numbers are implemented as Two's
 * complement and contains static_asserts failing if this is not the case.
 */
template <typename Int>
class integer_rank
{
    static_assert(std::is_integral<Int>::value,
                  "SignedInt has to be an integral type");

public:
    using int_type = Int;
    using rank_type = typename std::make_unsigned<int_type>::type;

    //! Maps value i to its rank in int_type. For any pair T x < y the invariant
    //! integer_rank<T>::rank_of_int(x) < integer_rank<T>::rank_of_int(y) holds.
    static constexpr rank_type rank_of_int(int_type i) {
        return use_identity_
               ? static_cast<rank_type>(i)
               : static_cast<rank_type>(i) ^ sign_bit_;
    }

    //! Returns the r-th smallest number of int_r. It is the inverse of
    //! rank_of_int, i.e. int_at_rank(rank_of_int(i)) == i for all i.
    static constexpr int_type int_at_rank(rank_type r) {
        return use_identity_
               ? static_cast<int_type>(r)
               : static_cast<int_type>(r ^ sign_bit_);
    }

private:
    constexpr static bool use_identity_ = !std::is_signed<int_type>::value;

    constexpr static rank_type sign_bit_
        = (rank_type(1) << (8 * sizeof(rank_type) - 1));

    // These test fail if a signed type does not use Two's complement
    static_assert(rank_of_int(std::numeric_limits<int_type>::min()) == 0,
                  "Rank of minimum is not zero");
    static_assert(rank_of_int(std::numeric_limits<int_type>::min() + 1) == 1,
                  "Rank of minimum+1 is not one");
    static_assert(rank_of_int(std::numeric_limits<int_type>::max())
                  == std::numeric_limits<rank_type>::max(),
                  "Rank of maximum is not maximum rank");
    static_assert(rank_of_int(std::numeric_limits<int_type>::max()) >
                  rank_of_int(int_type(0)),
                  "Rank of maximum is not larger than rank of zero");
};

//! Internal implementation of bitarray; do not invoke directly
//! \tparam Size  Number of bits the data structure is supposed to store
//! \tparam SizeIsAtmost64  Switch between inner node implementation (false)
//!                         and leaf implementation (true)
template <size_t Size, bool SizeIsAtmost64>
class bitarray_recursive;

template <size_t Size>
class bitarray_recursive<Size, false>
{
    static constexpr size_t leaf_width = 6;
    static constexpr size_t width = tlx::Log2<Size>::ceil;
    static_assert(width > leaf_width,
                  "Size has to be larger than 2**leaf_width");
    static constexpr size_t root_width = (width % leaf_width)
                                         ? (width % leaf_width)
                                         : leaf_width;
    static constexpr size_t child_width = width - root_width;
    using child_type = bitarray_recursive<1llu << child_width, child_width <= 6>;

    static constexpr size_t root_size = (Size + child_type::size - 1) / child_type::size; // div_ceil
    using root_type = bitarray_recursive<root_size <= 32 ? 32 : 64, true>;

    using child_array_type = std::array<child_type, root_size>;

public:
    static constexpr size_t size = Size;

    explicit bitarray_recursive() noexcept = default;
    bitarray_recursive(const bitarray_recursive&) noexcept = default;
    bitarray_recursive(bitarray_recursive&&) noexcept = default;
    bitarray_recursive& operator = (const bitarray_recursive&) noexcept = default;
    bitarray_recursive& operator = (bitarray_recursive&&) noexcept = default;

    void set_bit(const size_t i) {
        const auto idx = get_index_(i);
        root_.set_bit(idx.first);
        children_[idx.first].set_bit(idx.second);
    }

    void clear_bit(const size_t i) {
        const auto idx = get_index_(i);
        children_[idx.first].clear_bit(idx.second);
        if (children_[idx.first].empty())
            root_.clear_bit(idx.first);
    }

    bool is_set(const size_t i) const {
        const auto idx = get_index_(i);
        return children_[idx.first].is_set(idx.second);
    }

    void clear_all() {
        root_.clear_all();
        for (auto& child : children_)
            child.clear_all();
    }

    bool empty() const {
        return root_.empty();
    }

    size_t find_lsb() const {
        assert(!empty());

        const size_t child_idx = root_.find_lsb();
        const size_t child_val = children_[child_idx].find_lsb();

        return child_idx * child_type::size + child_val;
    }

private:
    child_array_type children_;
    root_type root_;

    std::pair<size_t, size_t> get_index_(size_t i) const {
        assert(i < size);
        return { i / child_type::size, i % child_type::size };
    }
};

template <size_t Size>
class bitarray_recursive<Size, true>
{
    static_assert(Size <= 64, "Support at most 64 bits");
    using uint_type = typename std::conditional<
        Size <= 32, uint32_t, uint64_t>::type;

public:
    static constexpr size_t size = Size;

    explicit bitarray_recursive() noexcept : flags_(0) { }
    bitarray_recursive(const bitarray_recursive&) noexcept = default;
    bitarray_recursive(bitarray_recursive&&) noexcept = default;
    bitarray_recursive& operator = (const bitarray_recursive&) noexcept = default;
    bitarray_recursive& operator = (bitarray_recursive&&) noexcept = default;

    void set_bit(const size_t i) {
        assert(i < size);
        flags_ |= uint_type(1) << i;
    }

    void clear_bit(const size_t i) {
        assert(i < size);
        flags_ &= ~(uint_type(1) << i);
    }

    bool is_set(const size_t i) const {
        assert(i < size);
        return (flags_ & (uint_type(1) << i)) != 0;
    }

    void clear_all() {
        flags_ = 0;
    }

    bool empty() const {
        return !flags_;
    }

    size_t find_lsb() const {
        assert(!empty());
        return tlx::ffs(flags_) - 1;
    }

private:
    uint_type flags_;
};

/*!
 * A bitarray of fixed size supporting reading, setting, and clearing
 * of individual bits. The data structure is optimized to find the bit with
 * smallest index that is set (find_lsb).
 *
 * The bitarray is implemented as a search tree with a fan-out of up to 64.
 * It is thus very flat, and all operations but with the exception of clear_all
 * have a complexity of O(log_64(Size)) which is << 10 for all practical purposes.
 */
template <size_t Size>
class bitarray
{
    using impl_type = bitarray_recursive<Size, Size <= 64>;

public:
    static constexpr size_t size = Size;

    explicit bitarray() noexcept = default;
    bitarray(const bitarray&) noexcept = default;
    bitarray(bitarray&&) noexcept = default;
    bitarray& operator = (const bitarray&) noexcept = default;
    bitarray& operator = (bitarray&&) noexcept = default;

    //! Set the i-th bit to true
    void set_bit(const size_t i) {
        impl_.set_bit(i);
    }

    //! Set the i-th bit to false
    void clear_bit(const size_t i) {
        impl_.clear_bit(i);
    }

    //! Returns value of the i-th
    bool is_set(const size_t i) const {
        return impl_.is_set(i);
    }

    //! Sets all bits to false
    void clear_all() {
        impl_.clear_all();
    }

    //! True if all bits are false
    bool empty() const {
        return impl_.empty();
    }

    //! Finds the bit with smallest index that is set
    //! \Warning: If empty() is true, the result is undefined
    size_t find_lsb() const {
        return impl_.find_lsb();
    }

private:
    impl_type impl_;
};

template <unsigned Radix, typename Int>
class bucket_computation
{
    static_assert(std::is_unsigned<Int>::value, "Require unsigned integer");
    static constexpr unsigned radix_bits = tlx::Log2<Radix>::floor;

public:
    //! Return bucket index key x belongs to given the current insertion limit
    size_t operator () (const Int x, const Int insertion_limit) const {
        constexpr Int mask = (1u << radix_bits) - 1;

        assert(x >= insertion_limit);

        const auto diff = x ^ insertion_limit;
        if (!diff) return 0;

        const auto diff_in_bit = (8 * sizeof(Int) - 1) - tlx::clz(diff);

        const auto row = diff_in_bit / radix_bits;
        const auto bucket_in_row = ((x >> (radix_bits * row)) & mask) - row;

        const auto result = row * Radix + bucket_in_row;

        return result;
    }

    //! Return smallest key possible in bucket idx assuming insertion_limit==0
    Int lower_bound(const size_t idx) const {
        assert(idx < num_buckets);

        if (idx < Radix)
            return static_cast<Int>(idx);

        const size_t row = (idx - 1) / (Radix - 1);
        const auto digit = static_cast<Int>(idx - row * (Radix - 1));

        return digit << radix_bits * row;
    }

    //! Return largest key possible in bucket idx assuming insertion_limit==0
    Int upper_bound(const size_t idx) const {
        assert(idx < num_buckets);

        if (idx == num_buckets - 1)
            return std::numeric_limits<Int>::max();

        return lower_bound(idx + 1) - 1;
    }

private:
    constexpr static size_t num_buckets_(size_t bits) {
        return (bits >= radix_bits)
               ? (Radix - 1) + num_buckets_(bits - radix_bits)
               : (1 << bits) - 1;
    }

public:
    //! Number of buckets required given Radix and the current data type Int
    static constexpr size_t num_buckets =
        num_buckets_(std::numeric_limits<Int>::digits) + 1;
};




} // namespace radixheap_detail

//! \addtogroup tlx_data_structures
//! \{

/**
 * This class implements a monotonic integer min priority queue, more specific
 * a multi-level radix heap.
 *
 * Here, monotonic refers to the fact that the heap maintains an insertion limit
 * and does not allow the insertion of keys smaller than this limit. The
 * frontier is increased to the current minimum when invoking the methods
 * top_key, top_data, pop and swap_top_bucket. To query the currently smallest
 * item without updating the insertion limit use peak_top_key.
 *
 * We implement a two level radix heap. Let k=sizeof(KeyType)*8 be the number of
 * bits in a key. In contrast to an ordinary radix heap which contains k buckets,
 * we maintain ceil(k/log2(Radix)) rows each containing Radix-many buckets.
 * This reduces the number of move operations when reorganizing the data structure.
 *
 * The implementation loosly follows the description of "An Experimental Study
 * of Priority Queues in External Memory" [Bregel et al.] and is also inspired
 * by https://github.com/iwiwi/radix-heap
 *
 * \tparam KeyType   Has to be an unsigned integer type
 * \tparam DataType  Type of data payload
 * \tparam Radix     A power of two <= 64.
 */
template <typename ValueType, typename KeyExtract, typename KeyType, unsigned Radix = 8>
class radixheap
{
    static_assert(radixheap_detail::tlx::Log2<Radix>::floor == radixheap_detail::tlx::Log2<Radix>::ceil,
                  "Radix has to be power of two");

    static constexpr bool debug = false;

public:
    using key_type = KeyType;
    using value_type = ValueType;
    using bucket_index_type = size_t;

    static constexpr unsigned radix = Radix;

protected:
    using encoder = radixheap_detail::integer_rank<key_type>;
    using ranked_key_type = typename encoder::rank_type;
    using bucket_map_type = radixheap_detail::bucket_computation<Radix, ranked_key_type>;

    static constexpr unsigned radix_bits = radixheap_detail::tlx::Log2<radix>::floor;
    static constexpr unsigned num_layers = (8 * sizeof(ranked_key_type) + radix_bits - 1) / radix_bits; // div_ceil
    static constexpr unsigned num_buckets = bucket_map_type::num_buckets;

public:
    using bucket_data_type = std::deque<value_type>;


    explicit radixheap(KeyExtract key_extract = KeyExtract{})
        : key_extract_(key_extract) {
        initialize_();
    }

    // Copy
    radixheap(const radixheap&) = default;
    radixheap& operator = (const radixheap&) = default;

    // Move
    radixheap(radixheap&&) = default;
    radixheap& operator = (radixheap&&) = default;

    bucket_index_type get_bucket(const value_type& value) const {
        return get_bucket_key(key_extract_(value));
    }

    bucket_index_type get_bucket_key(const key_type key) const {
        const auto enc = encoder::rank_of_int(key);
        assert(enc >= insertion_limit_);

        return bucket_map_(enc, insertion_limit_);
    }

    //! Construct and insert element with priority key
    //! \warning In contrast to all other methods the key has to be provided
    //! explicitly as the first argument. All other arguments are passed to
    //! the constructor of the element.
    template <typename... Args>
    bucket_index_type emplace(const key_type key, Args&& ... args) {
        const auto enc = encoder::rank_of_int(key);
        assert(enc >= insertion_limit_);
        const auto idx = bucket_map_(enc, insertion_limit_);

        emplace_in_bucket(idx, std::forward<Args>(args)...);
        return idx;
    }

    //! Construct and insert element into bucket idx (useful if an item
    //! was inserted into the same bucket directly before)
    //! \warning Calling any method which updates the current
    //! can invalidate this hint
    template <typename... Args>
    void emplace_in_bucket(const bucket_index_type idx, Args&& ... args) {
        if (buckets_data_[idx].empty()) filled_.set_bit(idx);
        buckets_data_[idx].emplace_back(std::forward<Args>(args) ...);

        const auto enc = encoder::rank_of_int(key_extract_(buckets_data_[idx].back()));
        if (mins_[idx] > enc) mins_[idx] = enc;
        assert(idx == bucket_map_(enc, insertion_limit_));

        size_++;
    }

    //! Insert element with priority key
    bucket_index_type push(const value_type& value) {
        const auto enc = encoder::rank_of_int(key_extract_(value));
        assert(enc >= insertion_limit_);

        const auto idx = bucket_map_(enc, insertion_limit_);

        push_to_bucket(idx, value);

        return idx;
    }

    //! Insert element into specific bucket (useful if an item
    //! was inserted into the same bucket directly before)
    //! \warning Calling any method which updates the current
    //! can invalidate this hint
    void push_to_bucket(const bucket_index_type idx, const value_type& value) {
        const auto enc = encoder::rank_of_int(key_extract_(value));

        assert(enc >= insertion_limit_);
        assert(idx == get_bucket(value));

        if (buckets_data_[idx].empty()) filled_.set_bit(idx);
        buckets_data_[idx].push_back(value);

        if (mins_[idx] > enc) mins_[idx] = enc;

        size_++;
    }

    //! Indicates whether size() == 0
    bool empty() const {
        return size() == 0;
    }

    //! Returns number of elements currently stored
    size_t size() const {
        return size_;
    }

    //! Returns currently smallest key without updating the insertion limit
    key_type peak_top_key() const {
        assert(!empty());
        const auto first = filled_.find_lsb();
        return encoder::int_at_rank(mins_[first]);
    }

    //! Returns currently smallest key and data
    //! \warning Updates insertion limit; no smaller keys can be inserted later
    const value_type& top() {
        reorganize_();
        return buckets_data_[current_bucket_].back();
    }

    //! Removes smallest element
    //! \warning Updates insertion limit; no smaller keys can be inserted later
    void pop() {
        reorganize_();
        buckets_data_[current_bucket_].pop_back();
        if (buckets_data_[current_bucket_].empty())
            filled_.clear_bit(current_bucket_);
        --size_;
    }

    //! Exchanges the top buckets with an *empty* user provided bucket.
    //! Can be used for bulk removals and may reduce allocation overhead
    //! \warning The exchange bucket has to be empty
    //! \warning Updates insertion limit; no smaller keys can be inserted later
    void swap_top_bucket(bucket_data_type& exchange_bucket) {
        reorganize_();

        assert(exchange_bucket.empty());
        size_ -= buckets_data_[current_bucket_].size();
        buckets_data_[current_bucket_].swap(exchange_bucket);

        filled_.clear_bit(current_bucket_);
    }

    //! Clears all internal queues and resets insertion limit
    void clear() {
        for (auto& x : buckets_data_) x.clear();
        initialize_();
    }

protected:
    KeyExtract key_extract_;
    size_t size_ { 0 };
    ranked_key_type insertion_limit_{ 0 };
    size_t current_bucket_{ 0 };

    bucket_map_type bucket_map_;

    std::array<bucket_data_type, num_buckets> buckets_data_;

    std::array<ranked_key_type, num_buckets> mins_;
    radixheap_detail::bitarray<num_buckets> filled_;

    void initialize_() {
        size_ = 0;
        insertion_limit_ = std::numeric_limits<ranked_key_type>::min();
        current_bucket_ = 0;

        std::fill(mins_.begin(), mins_.end(), std::numeric_limits<ranked_key_type>::max());

        filled_.clear_all();
    }

    void reorganize_() {
        assert(!empty());

        // nothing do to if we already know a suited bucket
        if ((!buckets_data_[current_bucket_].empty())) {
            assert(current_bucket_ < Radix);
            return;
        }

        // mark current bucket as empty
        mins_[current_bucket_] = std::numeric_limits<ranked_key_type>::max();
        filled_.clear_bit(current_bucket_);

        // find a non-empty bucket
        const auto first_non_empty = filled_.find_lsb();
        #ifndef NDEBUG
        {
            assert(first_non_empty < num_buckets);

            for (size_t i = 0; i < first_non_empty; i++) {
                assert(buckets_data_[i].empty());
                assert(mins_[i] == std::numeric_limits<ranked_key_type>::max());
            }

            assert(!buckets_data_[first_non_empty].empty());
        }
        #endif

        if ((first_non_empty < Radix)) {
            // the first_non_empty non-empty bucket belongs to the smallest row
            // it hence contains only one key and we do not need to reorganise
            current_bucket_ = first_non_empty;
            return;
        }

        // update insertion limit
        {
            const auto new_ins_limit = mins_[first_non_empty];
            assert(new_ins_limit > insertion_limit_);
            insertion_limit_ = new_ins_limit;
        }

        auto& data_source = buckets_data_[first_non_empty];

        for(; !data_source.empty(); data_source.pop_front()) {
            auto& x = data_source.front();

            const ranked_key_type key = encoder::rank_of_int(key_extract_(x));
            assert(key >= mins_[first_non_empty]);
            assert(first_non_empty == mins_.size() - 1 || key < mins_[first_non_empty + 1]);
            const auto idx = bucket_map_(key, insertion_limit_);
            assert(idx < first_non_empty);

            // insert into bucket
            if (buckets_data_[idx].empty()) filled_.set_bit(idx);
            buckets_data_[idx].push_back(std::move(x));
            if (mins_[idx] > key) mins_[idx] = key;
        }

        data_source.clear();

        // mark consumed bucket as empty
        mins_[first_non_empty] = std::numeric_limits<ranked_key_type>::max();
        filled_.clear_bit(first_non_empty);

        // update global pointers and minima
        current_bucket_ = filled_.find_lsb();
        assert(current_bucket_ < Radix);
        assert(!buckets_data_[current_bucket_].empty());
        assert(mins_[current_bucket_] >= insertion_limit_);
    }
};

/**
 * Helper to easily derive type of radixheap for a pre-C++17 compiler.
 * Refer to radixheap for description of parameters.
 * \note Expects DataType to be default constructable
 */
template <typename DataType, unsigned Radix = 8, typename KeyExtract = void>
auto make_radixheap(KeyExtract&& key_extract)->radixheap < DataType, KeyExtract, decltype(key_extract(DataType { })), Radix >
{
    return (radixheap < DataType, KeyExtract, decltype(key_extract(DataType{ })), Radix > {
                key_extract
            });
}

//! \}

} // namespace Aux

#endif // !TLX_CONTAINER_RADIXHEAP_HEADER

/******************************************************************************/
