from libcpp.utility cimport pair

cdef extern from "<unordered_map>" namespace "std":
    cdef cppclass unordered_map[T, U]:
        cppclass iterator:
            pair[T, U]& operator*() nogil
            iterator operator++() nogil
            iterator operator--() nogil
            bint operator==(iterator) nogil
            bint operator!=(iterator) nogil
        #cppclass const_iterator(iterator):
        #    pass
        unordered_map() nogil except +
        unordered_map(unordered_map&) nogil except +
        #unordered_map(key_compare&)
        U& operator[](T&) nogil
        #unordered_map& operator=(unordered_map&)
        bint operator==(unordered_map&, unordered_map&) nogil
        bint operator!=(unordered_map&, unordered_map&) nogil
        bint operator<(unordered_map&, unordered_map&) nogil
        bint operator>(unordered_map&, unordered_map&) nogil
        bint operator<=(unordered_map&, unordered_map&) nogil
        bint operator>=(unordered_map&, unordered_map&) nogil
        U& at(T&) nogil
        iterator begin() nogil
        #const_iterator begin()
        size_t bucket(T&)
        size_t bucket_count()
        size_t bucket_size(size_t)
        void clear() nogil
        size_t count(T&) nogil
        #emplace?
        bint empty() nogil
        iterator end() nogil
        #const_iterator end()
        pair[iterator, iterator] equal_range(T&) nogil
        #pair[const_iterator, const_iterator] equal_range(key_type&)
        iterator erase(iterator) nogil
        iterator erase(iterator, iterator) nogil
        size_t erase(T&) nogil
        iterator find(T&) nogil
        #const_iterator find(key_type&)
        #get_allocator?
        #hasher hash_function?
        pair[iterator, bint] insert(U&)
        iterator insert(iterator,U&)
        pair[iterator, bint] insert(pair[T, U]) nogil # XXX pair[T,U]&  #single element
        iterator insert(iterator, pair[T, U]) nogil # XXX pair[T,U]&	#with hint
        #void insert(input_iterator, input_iterator) maybe this would work?
    	#insert(P&&) ?
        #key_equal key_eq()
        float load_factor()
        size_t max_bucket_count()
        float max_load_factor()
        void max_load_factor(float)
        size_t max_size() nogil
        void rehash(size_t)
        void reserve(size_t)
        size_t size() nogil
        void swap(unordered_map&) nogil
