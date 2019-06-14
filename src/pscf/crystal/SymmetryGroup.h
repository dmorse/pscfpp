#ifndef PSCF_SYMMETRY_GROUP_H
#define PSCF_SYMMETRY_GROUP_H

/*
* Pscfatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <vector>

namespace Pscf 
{

   using namespace Util;

   /**
   * Class template for a group of elements of type Symmetry.
   *
   * This is written as a template to allow the creation of groups that
   * use different types of objects to represent symmetry elements. The
   * simplest distinction is between point groups and full space groups.
   *
   * The algorithm requires only the template parameter class Symmetry 
   * satisfy the following requirements:
   *
   *  1) A Symmetry must be default constructible.
   *  2) An operator * is provided to represent element multiplication.
   *  3) Operators == and != are provided to represent equality & inequality.
   *  4) A method Symmetry::inverse() must return the inverse of a Symmetry. 
   *  5) A static method Symmetry::identity() must return the identity.
   *
   * \ingroup Crystal_Module
   */
   template <class Symmetry>    
   class SymmetryGroup 
   {

   public:

      /**
      * Default constructor.
      */
      SymmetryGroup();

      /**
      * Copy constructor.
      */
      SymmetryGroup(const SymmetryGroup& other);

      /**
      * Destructor.
      */
      ~SymmetryGroup();

      /**
      * Add a new element to the group.
      *
      * Return false if the element was already present, true otherwise.
      *
      * \param symmetry new symmetry element.
      * \return true if this is a new element, false if already present.
      */
      bool add(Symmetry& symmetry);

      /**
      * Generate a complete group from the current elements.
      */
      void makeCompleteGroup();

      /**
      * Find a symmetry within a group.
      *
      * Return a pointer to a symmetry if it is in the group,
      * or a null pointer if it is not.
      */
      const Symmetry* find(const Symmetry& symmetry) const;

      /**
      * Return a reference to the identity element.
      */
      const Symmetry& identity() const;

      /**
      * Return number of elements in group.
      */
      int size() const;

      /**
      * Assignment operator.
      */
      SymmetryGroup& operator = (const SymmetryGroup& other);
 
      /**
      * Element access operator (by reference).
      */
      Symmetry& operator [] (int i);
 
      /**
      * Element access operator (by reference).
      */
      const Symmetry& operator [] (int i) const;

      /**
      * Return true if valid complete group, or throw an Exception.
      */
      bool isValid() const;
 
   private:

      std::vector<Symmetry> elements_;

      Symmetry identity_;

   };


   // Method definitions

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
   const Symmetry* SymmetryGroup<Symmetry>::find(const Symmetry& symmetry) const
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
      const Symmetry* ptr = find(symmetry);
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

   /*
   * Return the current size of the group.
   */
   template <class Symmetry>
   inline int SymmetryGroup<Symmetry>::size() const
   {  return elements_.size(); }
 
   /*
   * Return identity element.
   */
   template <class Symmetry>
   inline const Symmetry& SymmetryGroup<Symmetry>::identity() const
   {  return identity_; }
 
   /*
   * Element access operator (by reference).
   */
   template <class Symmetry> inline 
   Symmetry& SymmetryGroup<Symmetry>::operator [] (int i)
   { return elements_[i]; }
 
   /*
   * Element access operator (by reference).
   */
   template <class Symmetry> inline 
   const Symmetry& SymmetryGroup<Symmetry>::operator [] (int i) const
   { return elements_[i]; }
 
}
#endif
