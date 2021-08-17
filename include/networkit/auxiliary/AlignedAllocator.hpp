
#ifndef NETWORKIT_AUXILIARY_ALIGNED_ALLOCATOR_HPP_
#define NETWORKIT_AUXILIARY_ALIGNED_ALLOCATOR_HPP_

#include <stdexcept>

/**
 * Allocator for aligned data.
 *
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */
template <typename T, size_t Alignment>
class AlignedAllocator {
public:
    // The following will be the same for virtually all allocators.
    typedef T *pointer;
    typedef const T *const_pointer;
    typedef T &reference;
    typedef const T &const_reference;
    typedef T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    T *address(T &r) const { return &r; }

    const T *address(const T &s) const { return &s; }

    size_t max_size() const {
        // The following has been carefully written to be independent of
        // the definition of size_t and to avoid signed/unsigned warnings.
        return (static_cast<size_t>(0) - static_cast<size_t>(1)) / sizeof(T);
    }

    // The following must be the same for all allocators.
    template <typename U>
    struct rebind {
        typedef AlignedAllocator<U, Alignment> other;
    };

    bool operator!=(const AlignedAllocator &other) const { return !(*this == other); }

    void construct(T *const p, const T &t) const {
        void *const pv = static_cast<void *>(p);

        new (pv) T(t);
    }

    void destroy(T *const p) const { p->~T(); }

    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const AlignedAllocator &) const { return true; }

    // Default constructor, copy constructor, rebinding constructor, and destructor.
    // Empty for stateless allocators.
    AlignedAllocator() {}

    AlignedAllocator(const AlignedAllocator &) {}

    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment> &) {}

    ~AlignedAllocator() {}

    // The following will be different for each allocator.
    T *allocate(size_t n) const {
        // The return value of allocate(0) is unspecified.
        // Mallocator returns NULL in order to avoid depending
        // on malloc(0)'s implementation-defined behavior
        // (the implementation can define malloc(0) to return NULL,
        // in which case the bad_alloc check below would fire).
        // All allocators can return NULL in this case.
        if (n == 0)
            return NULL;

        // All allocators should contain an integer overflow check.
        // The Standardization Committee recommends that std::length_error
        // be thrown in the case of integer overflow.
        if (n > max_size())
            throw std::length_error("AlignedAllocator<T>::allocate() - Integer overflow.");

        // Mallocator wraps malloc().
        void *const pv = _mm_malloc(n * sizeof(T), Alignment);

        // Allocators should throw std::bad_alloc in the case of memory allocation failure.
        if (pv == NULL)
            throw std::bad_alloc();

        return static_cast<T *>(pv);
    }

    void deallocate(T *const p, size_t) const { _mm_free(p); }

    // The following will be the same for all allocators that ignore hints.
    template <typename U>
    T *allocate(size_t n, const U * /* const hint */) const {
        return allocate(n);
    }

    // Allocators are not required to be assignable, so
    // all allocators should have a private unimplemented
    // assignment operator. Note that this will trigger the
    // off-by-default (enabled under /Wall) warning C4626
    // "assignment operator could not be generated because a
    // base class assignment operator is inaccessible" within
    // the STL headers, but that warning is useless.
private:
    AlignedAllocator &operator=(const AlignedAllocator &);
};

#endif // NETWORKIT_AUXILIARY_ALIGNED_ALLOCATOR_HPP_
