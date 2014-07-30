from libcpp.utility cimport pair

cdef extern from "<unordered_set>" namespace "std":
    cdef cppclass unordered_set[T]:
        cppclass iterator:
            T& operator*()
            iterator operator++() nogil
            iterator operator--() nogil
            bint operator==(iterator) nogil
            bint operator!=(iterator) nogil
        #cppclass const_iterator(iterator):
        #    pass
        unordered_set() nogil except +
        unordered_set(unordered_set&) nogil except +
        #set(key_compare&)
        #set& operator=(set&)
        bint operator==(unordered_set&, unordered_set&) nogil
        bint operator!=(unordered_set&, unordered_set&) nogil
        bint operator<(unordered_set&, unordered_set&) nogil
        bint operator>(unordered_set&, unordered_set&) nogil
        bint operator<=(unordered_set&, unordered_set&) nogil
        bint operator>=(unordered_set&, unordered_set&) nogil
        iterator begin() nogil
        #const_iterator begin()
        size_t bucket(T&)
        size_t bucket_count() 
        size_t bucket_size(size_t n)
        void clear() nogil
        size_t count(T&) nogil
        #emplace?
        bint empty() nogil
        iterator end() nogil
        #const_iterator end()
        pair[iterator, iterator] equal_range(T&) nogil
        #pair[const_iterator, const_iterator] equal_range(T&)
        iterator erase(iterator) nogil
        iterator erase(iterator, iterator) nogil
        size_t erase(T&) nogil
        iterator find(T&) nogil
        #const_iterator find(T&)
        #get_allocator()?
        #hasher hash_function()?
        pair[iterator, bint] insert(T&) nogil
        iterator insert(iterator, T&) nogil
        #void insert(input_iterator, input_iterator)
        #key_eq key_equal()
        float load_factor()
        size_t max_bucket_count()
        float max_load_factor()
        void max_load_factor(float)
        size_t max_size() nogil
        void rehash(size_t)
        void reserve(size_t)
        size_t size() nogil
        void swap(unordered_set&) nogil
