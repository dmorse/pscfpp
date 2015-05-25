namespace Util
{

   /**
   * \defgroup Container_Module Container Templates
   * \ingroup Util_NS_Module
   *
   * A set of container and iterator class templates.
   *
   * This module contains a set of simple container templates, some of which
   * are similar to containers provided by the C++ standard library. Bounds 
   * checking of indices for all array containers can be turned on (for 
   * safety) or off (for speed) by defining or not defining the UTIL_DEBUG 
   * preprocessor macro.
   *
   * \section container_array_matrix_sec Array and Matrix Containers
   * 
   * Containers templates whose name contains the string 'Array' are one 
   * dimensional array containers, much like C arrays. All such containers
   * overload the subscript [] operator so as to return an object by 
   * reference, using the same syntax as a C array or a std::vector: If
   * A is an array, then A[i] is a reference to the ith element of A.
   * 
   * Container templates whose name contains the string 'Matrix' are two 
   * dimensional arrays. These overload the (int, int) operator to access
   * elements: If M is a Matrix, then M(i, j) is a reference to the element 
   * in column j of row i of A. 
   *
   * \section container_prefix_sec Container Name Prefixes
   *
   * The names of many containers have prefixes before the word Array or
   * Matrix that indicates policies for memory allocation and management.
   *
   * Containers templates whose name begins with the letter 'D' (such as
   * DArray, DSArray, DPArray, and DMatrix) use dynamically allocated memory. 
   * The declaration "DArray<int> A" declares a dynamically allocated array
   * of integers.  Memory must be explicitly allocated for these containers 
   * by calling the "allocate" method after the container is instantiated
   * and before it is used. Dynamically allocated containers can only be 
   * allocated once and are not resizable. Attempting to allocate a 
   * container more than once is as an error, and causes an Exception to
   * be thrown.
   * 
   * Containers templates whose name begins with a letter 'F' (such as
   * FArray, FSArray, FPArray, and FMatrix) are fixed size containers.
   * The capacity of each such container is determined at compile time 
   * by a template parameter or parameters. Thus, for example, 
   * \code
   *    FArray<int, 4> A;
   * \endcode
   * declares a fixed size array of four integers, much like the 
   * declaration "int V[4]" of a fixed size C array.
   *
   * The letter "S" in the names of DSArray and FSArray indicate that these
   * are "sized" arrays. These arrays have a variable logical size that is 
   * less than or equal to the physical capacity. The logical size is the 
   * current number of elements, which are always stored contiguously from
   * index 0 to index size - 1. Accessing an element with index greater than
   * or equal to size is an error, and will cause an Exception to be thrown
   * if debugging is enabled. The capacity of an array is the number of 
   * elements for which memory has been allocated. The size of such an array
   * is initially set to zero, and elements are added sequentially by the 
   * append() method, which adds a new element at the end of the array and 
   * increments the size counter. Once the size reaches the array capacity, 
   * attempting to append another element will cause an Exception to be 
   * thrown.
   *
   * Array containers whose name includes the prefix G are sized arrays
   * with a capacity that can grow (G="growable") as needed as elements 
   * are appended. The "GArray" template thus implements a dynamic array 
   * very similiar to the standard library std::vector. Automatic 
   * resizing changes the address of the beginning of the array, and
   * invalidates all iterators and pointers to elements. 
   *
   * \section container_pointer_sec Pointer Arrays
   *
   * Container templates whose name contains the prefix "P" are pointer
   * arrays. A pointer array is a container that stores pointers to objects 
   * that are instantiated outside of the array, rather than storing actual 
   * objects.  The containers DPArray and FPArray are dynamically allocated 
   * fixed size pointer arrays, respectively.  The GPArray array is a growable
   * pointer array, which can grow without bound. The pointer arrays are 
   * all similar to "sized" arrays in that they have a logical size that 
   * must be less than or equal to their capacity, and in that elements 
   * must can be added to the end of an initially empty array by function
   * "append(T& )". Pointer arrays use the same interface for the 
   * subscript operator (which returns a reference) and the append function 
   * (which takes a reference parameter) as that used by the sized and 
   * object arrays.  A pointer array of type DPArray< T > is thus different 
   * from a sized array of pointers, of type DSArray<T*>, because DPArray< T > 
   * overloads the [] operator to return a reference to an object of type T 
   * that is referenced by a private pointer, whereas the subscript operator
   * for a DSArray<T*> returns an actual pointer. 
   */

   /**
   * \defgroup Array_Module Object Arrays
   *
   * The Array containers that do not have a P prefix are one-dimensional 
   * array containers that store objects by value. These all overload the
   * subscript [] operator to provide access to elements as references.
   *
   * The DArray and FArray containers are simple wrappers for dynamically
   * allocated and fixed-size C arrays, respectively.  The DSArray and 
   * FSArray containers are dynamically and statically allocated arrays, 
   * respectively, that have both a fixed capacity but a variable logical 
   * size, with contiguous elements. The GArray container is a growable
   * sized array, similar to a std::vector. Destructors for these arrays
   * all delete the associated C array of objects.
   *
   * An RArray < T > is an Array that is intended to be used as an alias
   * for, or a shallow copy of, a target DArray, FArray or C array. An RArray
   * contains a copy of the array address and capacity of the target array,
   * where are copied by the RArray::associate() method. Like other array
   * containers, an RArray overloads the [] operator to provide access to 
   * elements as references. The destructor of an RArray does not delete 
   * the associated C array.
   *
   * A RingBuffer is a cylic buffer array for which the append() method
   * adds elements to the end of a sequence if the buffer is not full, and
   * overwrites the oldest element if it is.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup Pointer_Array_Module Pointer Arrays
   *
   * One-dimensional array containers that store pointers, and associated 
   * iterators.
   *
   * The DPArray and FPArray class templates are dynamically and statically
   * allocated pointer arrays, respectively. A GPArray is a growable pointer
   * array. 
   * 
   * Each DPArray < T >, FArray <T, N>, or GPArray<T> container has a private
   * C array of T* pointers. These containers all overload the the [] operator 
   * so as to return a T& reference, rather than a T* pointer. The append
   * method takes a T& reference as a parameter. The destructor for a pointer 
   * array deletes the underlying array of T* pointers, but not the T objects 
   * to which they point.
   * 
   * An ArrayStack < T > container is a finite capacity stack that is
   * implemented as a dynamically allocated array of T* pointers. Objects
   * can be pushed onto or popped off the top of the stack using the 
   * push(T&) and pop() methods. An ArrayStack can be allocated only once, 
   * and cannot be resized.
   *
   * An SSet< T > is a container that holds pointers to an unordered
   * set of T objects. It provides fast addition and removal of single
   * elements, and fast iteration through all elements. The indexing of 
   * elements is arbitrary and mutable, and may change when an element 
   * is deleted.
   *
   * An ArraySet < T > is a container that holds pointers to a subset of
   * the elements of an associated array of T objects. The indexing of
   * elements within an ArraySet container is arbitrary and mutable.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup Matrix_Module Matrix Containers
   *
   * Two-dimensional array containers that store by objects value.
   *
   * Matrix containers overload the () operator to return elements
   * by refernce. If A is a matrix object A(i, j) returns element
   * i, j of the matrix.
   *
   * The Matrix base class defines a conventional matrix, in which
   * all rows are the same length. The DMatrix and FMatrix subclasses
   * use dynamically allocated and fixed-size arrays, respectively.
   *
   * The RaggedMatrix base class and DRaggedMatrix subclass define
   * two-dimensional containers in which different rows can have
   * different lengths, though the list of row dimensions can be 
   * specified only once.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup List_Module Linked List
   *
   * A simple linked list implementation and associated iterator.
   *
   * \ingroup Container_Module
   */

   /**
   * \defgroup Iterator_Module Iterators
   *
   * Iterators for use with particular containers.
   *
   * \ingroup Container_Module
   */

}
