#ifndef PRDC_SYMMETRY_GROUP_TPP
#define PRDC_SYMMETRY_GROUP_TPP

/*
* Pscfatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SymmetryGroup.h"
#include <util/global.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Member function definitions (non-inline)

   /*
   * Default constructor
   */
   template <class Symmetry>
   SymmetryGroup<Symmetry>::SymmetryGroup()
   {
      identity_ = Symmetry::identity();
      elements_.push_back(identity_);
   }

   /*
   * Copy constructor
   */
   template <class Symmetry>
   SymmetryGroup<Symmetry>::SymmetryGroup(const SymmetryGroup& other)
   {
      elements_.clear();
      identity_ = other.identity();
      for (int i = 0; i < other.size(); ++i) {
         elements_.push_back(other.elements_[i]);
      }
   }

   /*
   * Destructor
   */
   template <class Symmetry>
   SymmetryGroup<Symmetry>::~SymmetryGroup()
   {}

   /*
   * Assignment operator.
   */
   template <class Symmetry>
   SymmetryGroup<Symmetry>& 
   SymmetryGroup<Symmetry>::operator = (const SymmetryGroup& other)
   {
      if (this != &other) {
         identity_ = other.identity();
         elements_.clear();
         for (int i = 0; i < other.size(); ++i) {
            elements_.push_back(other.elements_[i]);
         }
      }
      return *this;
   }

   /*
   * Find an element in the group, return const pointer to its address.
   *
   * Return a null pointer if symmetry is not a member of the group.
   */
   template <class Symmetry>
   Symmetry const * SymmetryGroup<Symmetry>::find(Symmetry const& symmetry) 
   const
   {
      for (int i=0; i < size(); ++i) {
         if (symmetry == elements_[i]) {
            return &(elements_[i]);
         }
      }
      // Return null pointer if not found
      return 0;
   }
 
   /*
   * Add a new element to the group.
   */
   template <class Symmetry>
   bool SymmetryGroup<Symmetry>::add(Symmetry& symmetry)
   {
      // Check if symmetry is already present
      const Symmetry* ptr = find(symmetry);

      // If not already present, add symmetry to this group
      bool added; 
      if (ptr == 0) {
         elements_.push_back(symmetry);
         added = true;
      } else {
         added = false;
      }
      return added;
   }

   /*
   * Create a complete group.
   */
   template <class Symmetry>
   void SymmetryGroup<Symmetry>::makeCompleteGroup()
   {
      Symmetry a, b, c;
      int      i, j, n;
      bool added, complete;

      // Add all inverses
      n = size();
      for (i = 0; i < n; ++i) {
         a = elements_[i].inverse();
         add(a);
      }

      // Add all products of existing elements, and their inverses
      complete = false;
      while (!complete) {
         complete = true;
         n = size();
         for (i = 0; i < n; ++i) {
            a = elements_[i];
            for (j = 0; j < n; ++j) {
               b = elements_[j];
               c = a*b;
               added = add(c);
               if (added) {
                  complete = false;
               }
               b = c.inverse(); 
               added = add(b);
               if (added) { 
                  complete = false;
               }
            } 
         } 
      } 

   }

   /*
   * Clear all elements except the identity.
   */
   template <class Symmetry>
   void SymmetryGroup<Symmetry>::clear()
   {
      elements_.clear();
      identity_ = Symmetry::identity();
      elements_.push_back(identity_);
   }

   /*
   * Determine if two space groups are equivalent.
   */
   template <class Symmetry>
   bool 
   SymmetryGroup<Symmetry>::operator == 
                              (SymmetryGroup<Symmetry> const & other) const
   {
      if (size() != other.size()) {
         return false; 
      } else {
         Symmetry const * ptr = 0;
         for (int i = 0; i < size(); ++i) {
            ptr = other.find(elements_[i]);
            if (ptr == 0) return false;
         }
         for (int i = 0; i < other.size(); ++i) {
            ptr = find(other.elements_[i]);
            if (ptr == 0) return false;
         }
      }
      return true;
   }

   /*
   * Determine if two space groups are inequivalent.
   */
   template <class Symmetry>
   bool 
   SymmetryGroup<Symmetry>::operator != 
                              (SymmetryGroup<Symmetry> const & other) const
   { return !(*this == other); }

   /*
   * Check validity of this group.
   */
   template <class Symmetry>
   bool SymmetryGroup<Symmetry>::isValid() const
   {
      Symmetry a, b, c;
      int      i, j, n;
 
      // Check for inclusion of identity element
      c = Symmetry::identity();
      if (find(c) == 0) {
         UTIL_THROW("Identity element is not in group");
      }

      // Check inverses, uniqueness, and completeness
      n = size();
      for (i = 0; i < n; ++i) {
         a = elements_[i];
         c = a.inverse();
         if (find(c) == 0) {
            UTIL_THROW("Inverse of element not in group");
         }
         for (j = 0; j < n; ++j) {
            b = elements_[j];
            if (i != j) {
               if (a == b) {
                  UTIL_THROW("An element of the group is not unique");
               }
            }
            c = a*b;
            if (find(c) == 0) {
               UTIL_THROW("Product of two element is not in group");
            }
         }
      }

      // If no Exceptions have been thrown, return true
      return true; 

   }

}
}
#endif
