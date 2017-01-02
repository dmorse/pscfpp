/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Grid.h"
#include "Dimension.h"
#include <util/global.h>

namespace Util
{

   Grid::Grid()
    : dimensions_(0),
      offsets_(0),
      size_(0)
   {
      IntVector dimensions(1);
      setDimensions(dimensions); 
   }

   Grid::Grid(const IntVector& dimensions)
    : dimensions_(0),
      offsets_(0),
      size_(0)
   {
      setDimensions(dimensions); 
   }

   void Grid::setDimensions(const IntVector& dimensions)
   { 
      int i;
      for (i = 0; i < Dimension; ++i) {
         if (dimensions[i] <= 0) {
            UTIL_THROW("Grid dimensions must be positive");
         }
      }
 
      dimensions_ = dimensions;
      offsets_[Dimension -1] = 1;
      for (i = Dimension - 1; i > 0; --i) {
         offsets_[i-1] = offsets_[i]*dimensions_[i];
      }
      size_ = offsets_[0]*dimensions_[0];
   }

   int Grid::rank(const IntVector& position) const
   {
      int i;
      int result = 0;
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

   IntVector Grid::position(int id) const
   {
      IntVector position;
      int       remainder = id;

      int i;
      for (i = 0; i < Dimension - 1; ++i) {
         position[i] = remainder/offsets_[i]; 
         remainder -= position[i]*offsets_[i]; 
      }
      position[i] = remainder;
      return position;
   }

   bool Grid::isInGrid(int coordinate, int i) const
   {
      bool result = true;
      if (coordinate <  0) 
         result = false;
      if (coordinate >= dimensions_[i]) 
         result = false;
      return result;
   }

   bool Grid::isInGrid(IntVector& position) const
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

   int Grid::shift(int& coordinate, int i) const
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

   IntVector Grid::shift(IntVector& position) const
   {
      IntVector shifts;
      for (int i = 0; i < Dimension; ++i) {
         shifts[i] = shift(position[i], i);
      }
      return shifts;
   }

} 
