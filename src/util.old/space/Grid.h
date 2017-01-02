#ifndef UTIL_GRID_H
#define UTIL_GRID_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/space/IntVector.h>

namespace Util
{

   /**
   * A grid of points indexed by integer coordinates.
   *
   * The coordinates of a point on a grid form an IntVector, referred 
   * to here as a grid position. Each element of a grid position must
   * lie in the range 0 <= position[i] < dimension(i), where i indexes
   * a Cartesian axis, and dimension(i) is the dimension of the grid 
   * along axis i. 
   * 
   * Each grid position is also assigned a non-negative integer rank.  
   * Grid position ranks are ordered sequentially like elements in 
   * a multi-dimensional C array, with the last coordinate being the 
   * most rapidly varying.
   *
   * \ingroup Space_Module
   */
   class Grid
   {

   public:

      /**
      * Default constructor
      */
      Grid();

      /**
      * Constructor
      *
      * \param dimensions  IntVector of grid dimensions
      */
      Grid(const IntVector& dimensions);

      // Compiler default destructor.

      // Compiler copy constructor.

      /**
      * Set the grid dimensions in all directions.
      *
      * \param dimensions  IntVector of grid dimensions.
      */
      void setDimensions(const IntVector& dimensions);

      /**
      * Get an IntVector of the grid dimensions.
      */
      IntVector dimensions() const;

      /**
      * Get grid dimension along Cartesian direction i.
      *
      * \param i  index of Cartesian direction 0 <=i < Util::Dimension
      */
      int dimension(int i) const;

      /**
      * Get total number of grid points.
      */
      int size() const;

      /**
      * Get the position IntVector of a grid point with a specified rank.
      *
      * \param rank  integer rank of a grid point.
      * \return IntVector containing coordinates of specified point.
      */
      IntVector position(int rank) const;

      /**
      * Get the rank of a grid point with specified position.
      *
      * \param position  integer position of a grid point
      * \return integer rank of specified grid point
      */
      int rank(const IntVector& position) const;

      /**
      * Is this coordinate in range?
      *
      * 
      * \param coordinate  coordinate value for direction i
      * \param i  index for Cartesian direction
      * \return true iff 0 <= coordinate < dimension(i).
      */
      bool isInGrid(int coordinate, int i) const;

      /**
      * Is this IntVector grid position within the grid?
      *
      * Returns true iff 0 <= coordinate[i] < dimension(i) for all i.
      *
      * \param position  grid point position
      * \return true iff 0 <= coordinate[i] < dimension(i) for all i.
      */
      bool isInGrid(IntVector& position) const;

      /**
      * Shift a periodic coordinate into range.
      *
      * Upon return, the coordinate will be shifted to lie within the
      * range 0 <= coordinate < dimension(i) by subtracting an integer
      * multiple of dimension(i), giving coordinate - shift*dimension(i). 
      * The return value is the required integer `shift'.
      *
      * \param coordinate  coordinate in Cartesian direction i.
      * \param i  index of Cartesian direction, i >= 0. 
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
      * \param position  IntVector position within a grid.
      * \return IntVector of integer shifts.
      */
      IntVector shift(IntVector& position) const;

   private:

      /// Dimensions of grid
      IntVector dimensions_;

      /// Number of elements per increment in each direction
      IntVector offsets_;

      /// Total number of grid points
      int size_;

   };

   // Inline member functions

   inline IntVector Grid::dimensions() const
   {  return dimensions_; }

   inline int Grid::dimension(int i) const
   {
      assert(i >=0);
      assert(i < Dimension);  
      return dimensions_[i]; 
   }

   inline int Grid::size() const
   {  return size_; }

} 
#endif
