# distutils: language=c++

cdef class Cover:
	""" 
	Cover(n=0)
	
	Implements a cover of a set, i.e. an assignment of its elements to possibly overlapping subsets.
	
	Parameters
	----------
	n : int or networkit.Partition, optional
		Used for initialization of the cover. Either a node or a partition. Default: 0
	"""
	def __cinit__(self, n=0):
		if isinstance(n, Partition):
			self._this = move(_Cover((<Partition>n)._this))
		else:
			self._this = move(_Cover(<count?>n))

	cdef setThis(self, _Cover& other):
		swap[_Cover](self._this, other)
		return self

	def subsetsOf(self, e):
		""" 
		subsetsOf(e)		

		Get the ids of subsets in which the element `e` is contained.

		Parameters
		----------
		e : int
			An element

		Returns
		-------
		list(int)
			A set of subset ids in which `e` is contained.
		"""
		return self._this.subsetsOf(e)

	def extend(self):
		"""
		extend()

		Add an additional element (node).

		Returns
		-------
		int
			Id of added node.
		"""
		return self._this.extend()

	def addToSubset(self, s, e):
		"""
		addToSubset(s, e)

		Add the (previously unassigned) element `e` to the set `s`.

		Parameters
		----------
		s : int
			The input subset.
		e : int
			The element to be added.
		"""
		self._this.addToSubset(s, e)

	def removeFromSubset(self, s, e):
		""" 
		removeFromSubset(s, e)
		
		Remove the element `e` from the set `s`.

		Parameters
		----------
		s : int
			The input subset.
		e : int
			The element to be removed.
		"""
		self._this.removeFromSubset(s, e)

	def moveToSubset(self, index s, index e):
		""" 
		moveToSubset(s, e)

		Move the element `e` to subset `s`, i.e. remove it from all other subsets and place it in the subset.

		Parameters
		----------
		s : int
			The input subset.
		e : int
			The element to be moved.
		"""
		self._this.moveToSubset(s, e)

	def toSingleton(self, index e):
		""" 
		toSingleton(e)

		Creates a singleton set containing the element `e` and returns the index of the new set.

		Parameters
		----------
		e : int
			The input element.

		Returns
		-------
		int
			The id of the new set.
		"""
		self._this.toSingleton(e)

	def allToSingletons(self):
		""" 
		allToSingletons()		

		Assigns every element to a singleton set. Set id is equal to element id. 
		"""
		self._this.allToSingletons()

	def mergeSubsets(self, index s, index t):
		""" 
		mergeSubsets(s, t)		

		Assigns the elements from both sets to a new set.

		Parameters
		----------
		s : int
			The first subset.
		t : int
			The second subset.
		"""
		self._this.mergeSubsets(s, t)

	def setUpperBound(self, index upper):
		"""
		setUpperBound(upper)

		Sets an upper bound for the subset ids that CAN be assigned.

		Parameters
		----------
		upper : int
			Upper bound.
		"""
		self._this.setUpperBound(upper)

	def upperBound(self):
		""" 
		upperBound()		

		Get an upper bound for the subset ids that have been assigned.
		(This is the maximum id + 1.)

		Returns
		-------
		int
			An upper bound.
		"""
		return self._this.upperBound()

	def lowerBound(self):
		""" 
		lowerBound()		

		Get a lower bound for the subset ids that have been assigned.

		Returns
		-------
		int
			A lower bound.
		"""
		return self._this.lowerBound()

	def contains(self, index e):
		"""
		contains(e)		

		Check if cover assigns a valid subset to the element `e`.

		Parameters
		----------
		e : int
			The input element.

		Returns
		-------
		bool
			True, if `e` is assigned to a valid subset, False otherwise.
		"""
		return self._this.contains(e)

	def inSameSubset(self, index e1, index e2):
		"""
		inSameSubset(e1, e2)
		
		Check if two elements `e1` and `e2` belong to the same subset.

		Parameters
		----------
	 	e1 : int
			The first element.
		e2 : int
			The second element.

		Returns
		-------
		bool
			True if `e1` and `e2` belong to the same subset; False otherwise.
		"""
		return self._this.inSameSubset(e1, e2)

	def subsetSizes(self):
		""" 
		subsetSizes()		

		Get a list of subset sizes.

		Returns
		-------
		list(int)
			A list of subset sizes.

		Notes
		-----
		Indices do not necessarily correspond to subset ids.
		"""
		return self._this.subsetSizes()

	def subsetSizeMap(self):
		""" 
		subsetSizeMap()		

		Get a map from subset id to size of the subset.

		Returns
		-------
		dict(int ``:`` int)
			A map from subset id to size of the subset.
		"""
		return self._this.subsetSizeMap()

	def getMembers(self, s):
		""" 
		getMembers(s)		

		Get the members of a specific subset `s`.

		Returns
		-------
		list(int)
			The list of members of subset `s`.
		"""
		return self._this.getMembers(s)

	def numberOfElements(self):
		""" 
		numberOfElements()		

		Get the current number of elements in this cover.

		Returns
		-------
		int
			The current number of elements.
		"""
		return self._this.numberOfElements()

	def numberOfSubsets(self):
		"""  
		numberOfSubsets()		

		Get the current number of sets in this cover.

		Returns
		-------
		int
			The number of sets in this cover.
		"""
		return self._this.numberOfSubsets()

	def getSubsetIds(self):
		""" 
		getSubsetIds()		

		Get the ids of nonempty subsets.

		Returns
		-------
		list(int)
			A list of ids of nonempty subsets.
		"""
		return self._this.getSubsetIds()

cdef class Partition:
	""" 
	Partition(z=0)

	Implements a partition of a set, i.e. a subdivision of the
	set into disjoint subsets.

	Create a new partition data structure for `z` elements.

	Parameters
	----------
	size : int, optional
		Maximum index of an element. Default: 0
	"""

	def __cinit__(self, index size=0, vector[index] data=[]):
		if data.size() != 0:
			self._this = move(_Partition(data))
		else:
			self._this = move(_Partition(size))

	def __len__(self):
		"""
		__len__()

		Returns
		-------
		int
			Number of elements in the partition.
		"""
		return self._this.numberOfElements()

	def __getitem__(self, index e):
		""" 
		__getitem__(index e)
		
		Get the set (id) in which the element `e` is contained.

		Parameters
		----------
	 	e : int
	 		Index of element.

		Returns
		-------
		int
			The index of the set in which `e` is contained.
		"""
		return self._this.subsetOf(e)

	def __setitem__(self, index e, index s):
		""" 
		__setitem__(e, s)

		Set the set (id) in which the element `e` is contained.

		Parameters
		----------
		e : int
			Index of the element
		s : int
			Index of the subset
		"""
		self._this.addToSubset(s, e)

	def __copy__(self):
		"""
		Generates a copy of the partition
		"""
		return Partition().setThis(_Partition(self._this))

	def __deepcopy__(self):
		"""
		Generates a copy of the partition
		"""
		return Partition().setThis(_Partition(self._this))

	cdef setThis(self,  _Partition& other):
		swap[_Partition](self._this,  other)
		return self

	def subsetOf(self, e):
		""" 
		subsetOf(e)
		
		Get the set (id) in which the element `e` is contained.

		Parameters
		----------
		e : int
			Index of element.

		Returns
		-------
		int
			The index of the set in which `e` is contained.
		"""
		return self._this.subsetOf(e)

	def extend(self):
		""" 
		extend()

		Extend the data structure and create a slot	for one more element.

		Initializes the entry to `none` and returns the index of the entry.

		Returns
		-------
		int
			The index of the new element.
		"""
		return self._this.extend()

	def addToSubset(self, s, e):
		"""
		addToSubset(s, e)		

		Add a (previously unassigned) element `e` to the set `s`.

		Parameters
		----------
		s : int
			The index of the subset.
		e : int
			The element to add.
		"""
		self._this.addToSubset(s, e)

	def moveToSubset(self, index s, index e):
		"""  
		moveToSubset(s, e)	

		Move the (previously assigned) element `e` to the set `s.

		Parameters
		----------
		s : int
			The index of the subset.
		e : int
			The element to move.
		"""
		self._this.moveToSubset(s, e)

	def toSingleton(self, index e):
		""" 
		toSingleton(e)		

		Creates a singleton set containing the element `e`.

		Parameters
		----------
		e : int
			The index of the element.
		"""
		self._this.toSingleton(e)

	def allToSingletons(self):
		""" 
		allToSingletons()		

		Assigns every element to a singleton set. Set id is equal to element id.
		"""
		self._this.allToSingletons()

	def mergeSubsets(self, index s, index t):
		""" 
		mergeSubsets(s, t)		

		Assigns the elements from both sets to a new set and returns the id of it.

		Parameters
		----------
		s : int
			Set to merge.
		t : int
			Set to merge.

		Returns
		-------
		int
			Id of newly created set.
		"""
		return self._this.mergeSubsets(s, t)

	def setUpperBound(self, index upper):
		"""
		setUpperBound(upper)		

		Sets an upper bound for the subset ids that **can** be assigned.

		Parameters
		----------
		upper : int
			Highest assigned subset id + 1.
		"""
		self._this.setUpperBound(upper)

	def upperBound(self):
		""" 
		upperBound()	

		Return an upper bound for the subset ids that have been assigned.
		(This is the maximum id + 1.)

		Returns
		-------
		int
			The upper bound.
		"""
		return self._this.upperBound()

	def lowerBound(self):
		""" 
		lowerBound()		

		Get a lower bound for the subset ids that have been assigned.

		Returns
		-------
		int
			The lower bound.
		"""
		return self._this.lowerBound()

	def compact(self, useTurbo = False):
		""" 
		compact(userTurbo=False)		

		Change subset IDs to be consecutive, starting at 0.

		Parameters
		----------
		useTurbo : bool, optional
			If set to True, the C++ core uses a vector instead of a map to assign new ids 
			which results in a shorter running time but possibly a large space overhead.
			Default: False
		"""
		self._this.compact(useTurbo)

	def contains(self, index e):
		""" 
		contains(e)		

		Check if partition assigns a valid subset to the element `e`.

		Parameters
		----------
		e : int
			The input element.

		Returns
		-------
		bool
			True if the assigned subset is valid; False otherwise.
		"""
		return self._this.contains(e)

	def inSameSubset(self, index e1, index e2):
		""" 
		inSameSubset(e1, e2)		

		Check if two elements `e1` and `e2` belong to the same subset.

		Parameters
		----------
		e1 : int
			The first Element.
		e2 : int
			The second Element.

		Returns
		-------
		bool
			True if `e1` and `e2` belong to same subset, False otherwise.
		"""
		return self._this.inSameSubset(e1, e2)

	def subsetSizes(self):
		""" 
		subsetSizes()		

		Get a list of subset sizes. Indices do not necessarily correspond to subset ids.

	 	Returns
	 	-------
	 	list(int)
	 		A list of subset sizes.
		"""
		return self._this.subsetSizes()

	def subsetSizeMap(self):
		""" 
		subsetSizeMap()		

		Get a map from subset id to size of the subset.

		Returns
		-------
		dict(int ``:`` int)
			A map from subset id to size of the subset.
		"""
		return self._this.subsetSizeMap()

	def getMembers(self, s):
		""" 
		getMembers(s)		

		Get the members of the subset `s`.

		Parameters
		----------
		s : int
			The input subset.

		Returns
		-------
		list(int)
			A list containing the members of `s`.
		"""
		return self._this.getMembers(s)

	def numberOfElements(self):
		"""
		numberOfElements()

		Returns
		-------
		int
			Number of elements in the partition.
		"""
		return self._this.numberOfElements()

	def numberOfSubsets(self):
		"""
		numberOfSubsets()

		Get the current number of sets in this partition.

		Returns
		-------
		int
			The current number of sets.
		"""
		return self._this.numberOfSubsets()

	def getVector(self):
		""" 
		getVector()		

		Get the actual vector representing the partition data structure.

		Returns
		-------
		list(int)
			List containing information about partitions.
		"""
		return self._this.getVector()

	def setName(self, string name):
		"""
		setName(name)

		Set a human-readable identifier `name` for the instance.

		Parameters
		----------
		name : str
			The input name.
		"""
		self._this.setName(name)

	def getName(self):
		"""
		getName()

		Get the human-readable identifier.

		Returns
		-------
		str
			The name of this partition.
		"""
		return self._this.getName()

	def getSubsetIds(self):
		""" 
		getSubsetIds()

		Get the ids of nonempty subsets.

		Returns
		-------
		list(int)
			A set of ids of nonempty subsets.
		"""
		return self._this.getSubsetIds()

	def __eq__(self, Partition other not None):
		"""
		__eq__(other)

		Compare self to other partition.

		Equality is independent of the used partition
		ids. This tries to construct a mapping between the
		partition ids and returns True if such a mapping can
		be constructed.

		Parameters
		----------
		other : networkit.Partition
			The partition to compare to.

		Returns
		-------
		bool
			If the partitions are equal.
		"""
		if self._this.numberOfElements() != other._this.numberOfElements():
			return False

		cdef index i = 0
		cdef dict selfToOther = dict()
		for index in range(self._this.numberOfElements()):
			selfSubset = self[i]
			if selfSubset in selfToOther:
				if selfToOther[selfSubset] != other[i]:
					return False
			else:
				selfToOther[selfSubset] = other[i]
		return True

