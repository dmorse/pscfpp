#ifndef UTIL_GRID_ARRAY_H
#define UTIL_GRID_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * Multi-dimensional array with the dimensionality of space.
   *
   * The memory for a GridArray is stored in a single one-dimensional C array.
   * The subscript [] operator is overloaded to return an element indexed by
   * a one-dimensional rank, and the () operator is overloaded to return an
   * element indexed by an IntVector of grid coordinates.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class GridArray
   {

   public:

      /**
      * Constructor.
      */
      GridArray();

      /**
      * Copy constructor.
      */
      GridArray(const GridArray<Data>& other);

      /**
      * Destructor.
      *
      * Delete dynamically allocated C array, if allocated.
      */
      ~GridArray();

      /**
      * Assignment.
      */
      GridArray<Data>& operator = (const GridArray<Data>& other);

      // Initialization

      /**
      * Allocate memory for a matrix.
      *
      * \param dimensions IntVector containing dimensions
      */
      void allocate(const IntVector& dimensions);

      /**
      * Serialize a GridArray to/from an Archive.
      *
      * \param ar  archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Return true if the GridArray has been allocated, false otherwise.
      */
      bool isAllocated() const;

      // Grid Interface

      /**
      * Get all dimensions of array as an IntVector.
      *
      * \return IntVector containing the number of elements in each direction.
      */
      const IntVector& dimensions();

      /**
      * Get number of grid points along direction i.
      *
      * \param i index of Cartesian direction 0 <=i < 3.
      */
      int dimension(int i) const;

      /**
      * Get total number of grid points.
      */
      int size() const;

      /**
      * Get the position IntVector of a grid point with a specified rank.
      *
      * \param  rank integer rank of a grid point.
      * \return IntVector containing coordinates of specified point.
      */
      IntVector position(int rank) const;

      /**
      * Get the rank of a grid point with specified position.
      *
      * \param  position integer position of a grid point
      * \return integer rank of specified grid point
      */
      int rank(const IntVector& position) const;

      /**
      * Is this 1D coordinate in range?
      *
      * Returns true iff 0 <= coordinate < dimension(i).
      * \param coordinate coordinate value for direction i
      * \param i          index for Cartesian direction
      */
      bool isInGrid(int coordinate, int i) const;

      /**
      * Is this position within the grid?
      *
      * Returns true iff 0 <= coordinate[i] < dimension(i) for all i.
      *
      * \param position grid point position
      */
      bool isInGrid(IntVector& position) const;

      /**
      * Shift a periodic 1D coordinate into primary range.
      *
      * Upon return, the coordinate will be shifted to lie within the
      * range 0 <= coordinate < dimension(i) by subtracting an integer
      * multiple of dimension(i), giving coordinate - shift*dimension(i).
      * The return value is the required integer `shift'.
      *
      * \param  coordinate  coordinate in Cartesian direction i.
      * \param  i           index of Cartesian direction, i >= 0.
      * \return multiple of dimension(i) subtracted from input value.
      */
      int shift(int& coordinate, int i) const;

      /**
      * Shift a periodic position into primary grid.
      *
      * Upon return, each element of the parameter position is shifted
      * to lie within the range 0 <= position[i] < dimension(i) by
      * adding or subtracting an integer multiple of dimension(i). The
      * IntVector of shift values is returned.
      *
      * \param  position IntVector position within a grid.
      * \return IntVector of integer shifts.
      */
      IntVector shift(IntVector& position) const;

      // Array Interface

      /**
      * Return element by const reference, indexed by 1D rank.
      *
      * \param rank 1D array index of element
      */
      const Data& operator[] (int rank) const;

      /**
      * Return element by reference, indexed by 1D rank.
      *
      * \param rank 1D rank of element
      */
      Data& operator[] (int rank);

      /**
      * Return element by const reference, indexed by IntVector position.
      *
      * \param position IntVector of coordinates.
      */
      const Data& operator() (const IntVector& position) const;

      /**
      * Return element by reference, indexed by IntVector position.
      *
      * \param position IntVector of coordinates.
      */
      Data& operator() (const IntVector& position);

   private:

      /// Pointer to 1D C array of all elements.
      Data*  data_;

      /// Number of elements per increment in each direction
      IntVector offsets_;

      /// Dimensions of grid
      IntVector dimensions_;

      /// Total number of grid points
      int size_;

   };

   // Method definitions

   /**
   * Constructor (protected).
   */
   template <typename Data>
   inline GridArray<Data>::GridArray()
    : data_(0),
      offsets_(),
      dimensions_(),
      size_(0)
   {}

   /*
   * Destructor.
   *
   * Delete dynamically allocated C array.
   */
   template <typename Data>
   GridArray<Data>::~GridArray()
   {
      if (data_) {
         Memory::deallocate<Data>(data_, size_);
         size_ = 0;
      }
   }

   /*
   * Copy constructor.
   */
   template <typename Data>
   GridArray<Data>::GridArray(const GridArray<Data>& other)
    : data_(0),
      offsets_(),
      dimensions_(),
      size_(0)
   {
      // Precondition
      if (other.data_ == 0) {
         UTIL_THROW("Other GridArray must be allocated");
      }
      if (isAllocated()) {
         UTIL_THROW("GridArray already allocated in copy constructor");
      }

      allocate(other.dimensions_);
      if (offsets_ != other.offsets_ ) {
         UTIL_THROW("Unequal offsets");
      }
      if (size_ != other.size_ ) {
         UTIL_THROW("Unequal sizes");
      }
      for (int i = 0; i < size_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Assignment.
   */
   template <typename Data>
   GridArray<Data>& GridArray<Data>::operator = (const GridArray<Data>& other)
   {
      // Check for self assignment.
      if (this == &other) {
         return *this;
      }

      // Precondition
      if (other.data_ == 0) {
         UTIL_THROW("Other GridArray must be allocated before assignment");
      }

      // If this GridArray if not allocated, allocate now.
      // If it is allocated, check that dimensions are equal.
      if (!isAllocated()) {
         allocate(other.dimensions_);
      } else {
         if (dimensions_ != other.dimensions_ ) {
            UTIL_THROW("Unequal dimensions");
         }
      }
      if (offsets_ != other.offsets_ ) {
         UTIL_THROW("Unequal offsets");
      }
      if (size_ != other.size_ ) {
         UTIL_THROW("Unequal sizes");
      }

      // Copy elements
      for (int i = 0; i < size_; ++i) {
         data_[i] = other.data_[i];
      }

      return *this;
   }

   /*
   * Set dimensions and allocate memory.
   */
   template <typename Data>
   void GridArray<Data>::allocate(const IntVector& dimensions)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a GridArray");
      }
      for (int i = 0; i < Dimension; ++i) {
         if (dimensions[i] <= 0) {
            UTIL_THROW("Dimension not positive");
         }
      }
      dimensions_ = dimensions;
      offsets_[Dimension -1] = 1;
      for (int i = Dimension - 1; i > 0; --i) {
         offsets_[i-1] = offsets_[i]*dimensions_[i];
      }
      size_ = offsets_[0]*dimensions_[0];
      Memory::allocate<Data>(data_, size_);
   }

   /*
   * Serialize a GridArray to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void GridArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      IntVector dimensions;
      if (Archive::is_saving()) {
         dimensions = dimensions_;
      }
      ar & dimensions;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            allocate(dimensions);
         }
      }
      for (int i = 0; i < size_; ++i) {
         ar & data_[i];
      }
   }

   /*
   * Get IntVector of dimensions.
   */
   template <typename Data>
   inline const IntVector& GridArray<Data>::dimensions()
   {  return dimensions_; }

   /*
   * Get dimension in direction i.
   */
   template <class Data>
   inline int GridArray<Data>::dimension(int i) const
   {  return dimensions_[i]; }

   /*
   * Get total number of grid points.
   */
   template <class Data>
   inline int GridArray<Data>::size() const
   {  return size_; }

   /*
   * Calculate 1D rank from grid coordinates.
   */
   template <typename Data>
   #ifdef UTIL_DEBUG
   int GridArray<Data>::rank(const IntVector& position) const
   {
      int result = 0;
      int i;
      for (i = 0; i < Dimension - 1; ++i) {
         assert(position[i] >= 0);
         assert(position[i] < dimensions_[i]);
         result += position[i]*offsets_[i];
      }
      assert(position[i] >= 0);
      assert(position[i] < dimensions_[i]);
      result += position[i];
      return result;
   }
   #else
   inline int GridArray<Data>::rank(const IntVector& position) const
   {
      return (position[0]*offsets_[0] + position[1]*offsets_[1] + position[2]);
   }
   #endif

   /*
   * Calculate coordinates from 1D rank.
   */
   template <typename Data>
   IntVector GridArray<Data>::position(int rank) const
   {
      IntVector position;
      int remainder = rank;

      int i;
      for (i = 0; i < Dimension - 1; ++i) {
         position[i] = remainder/offsets_[i];
         remainder -= position[i]*offsets_[i];
      }
      position[i] = remainder;
      return position;
   }

   /*
   * Test if a single coordinate is within range of grid.
   */
   template <typename Data>
   bool GridArray<Data>::isInGrid(int coordinate, int i) const
   {
      bool result = true;
      if (coordinate <  0)
         result = false;
      if (coordinate >= dimensions_[i])
         result = false;
      return result;
   }

   /*
   * Test if a IntVector position is within primary grid domain.
   */
   template <typename Data>
   bool GridArray<Data>::isInGrid(IntVector& position) const
   {
      bool result = true;
      for (int i = 0; i < Dimension; ++i) {
         if (position[i] <  0)
            result = false;
         if (position[i] >= dimensions_[i])
            result = false;
      }
      return result;
   }

   /*
   * Shift a 1D coordinate to primary domain.
   */
   template <typename Data>
   int GridArray<Data>::shift(int& coordinate, int i) const
   {
      int shift;
      if (coordinate >= 0) {
         shift = coordinate/dimensions_[i];
      } else {
         shift = -1 + ((coordinate+1)/dimensions_[i]);
      }
      coordinate -= shift*dimensions_[i];
      return shift;
   }

   /*
   * Shift a IntVector position to primary domain.
   */
   template <typename Data>
   IntVector GridArray<Data>::shift(IntVector& position) const
   {
      IntVector shifts;
      for (int i = 0; i < Dimension; ++i) {
         shifts[i] = shift(position[i], i);
      }
      return shifts;
   }

   /*
   * Return element by const reference, indexed by 1D rank.
   */
   template <typename Data>
   inline const Data& GridArray<Data>::operator[] (int rank) const
   {  return *(data_ + rank); }

   /*
   * Return element by reference, indexed by 1D rank.
   */
   template <typename Data>
   inline Data& GridArray<Data>::operator[] (int rank)
   {  return *(data_ + rank); }

   /*
   * Return element by const reference, indexed by IntVector of coordinates
   */
   template <typename Data>
   inline
   const Data& GridArray<Data>::operator() (const IntVector& position) const
   {  return *(data_ + rank(position)); }

   /*
   * Return element by reference, indexed by IntVector of coordinates
   */
   template <typename Data>
   inline Data& GridArray<Data>::operator() (const IntVector& position)
   {  return *(data_ + rank(position)); }

   /*
   * Return true if the GridArray has been allocated, false otherwise.
   */
   template <class Data>
   inline bool GridArray<Data>::isAllocated() const
   {  return (bool)(data_ != 0); }

}
#endif
