#ifndef PRDC_SYMMETRY_GROUP_H
#define PRDC_SYMMETRY_GROUP_H

/*
* Pscfatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

//#include <iostream>
#include <vector>

namespace Pscf {
namespace Prdc {

   /**
   * Class template for a group of elements.
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
   * \ingroup Prdc_Crystal_Module
   */
   template <class Symmetry>    
   class SymmetryGroup 
   {

   public:

      /**
      * Default constructor.
      *
      * After construction, the group contains only the identity element.
      */
      SymmetryGroup();

      /**
      * Copy constructor.
      */
      SymmetryGroup(const SymmetryGroup<Symmetry>& other);

      /**
      * Destructor.
      */
      ~SymmetryGroup();

      /**
      * Assignment operator.
      */
      SymmetryGroup<Symmetry>& 
      operator = (const SymmetryGroup<Symmetry>& other);
 
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
      * Remove all elements except the identity.
      *
      * Return group to its state after default construction.
      */
      void clear();

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
      * Return number of elements in group (i.e., the order of the group).
      */
      int size() const;

      /**
      * Element access operator (by reference).
      *
      * \param i  integer id for a symmetry element.
      */
      Symmetry& operator [] (int i);
 
      /**
      * Element access operator (by reference).
      *
      * \param i  integer id for a symmetry element.
      */
      const Symmetry& operator [] (int i) const;

      /**
      * Equality operator.
      *
      * Return true if this and other are equivalent symmetry groups,
      * false otherwise.
      *
      * \param  other the group to which this one is compared.
      */
      bool operator == (SymmetryGroup<Symmetry> const &  other) const;

      /**
      * Inequality operator.
      *
      * Return true if this and other are inequivalent symmetry groups,
      * and false if they are equivalent.
      *
      * \param  other the group to which this one is compared.
      */
      bool operator != (SymmetryGroup<Symmetry> const &  other) const;

      /**
      * Return true if valid complete group, or throw an Exception.
      */
      bool isValid() const;
 
   private:

      // Array of symmetry elements
      std::vector<Symmetry> elements_;

      // Identity element
      Symmetry identity_;

   };


   // Inline member function definitions

   /*
   * Return the current size of the group.
   */
   template <class Symmetry>
   inline 
   int SymmetryGroup<Symmetry>::size() const
   {  return elements_.size(); }
 
   /*
   * Return identity element.
   */
   template <class Symmetry>
   inline 
   const Symmetry& SymmetryGroup<Symmetry>::identity() const
   {  return identity_; }
 
   /*
   * Element access operator (by reference).
   */
   template <class Symmetry> 
   inline 
   Symmetry& SymmetryGroup<Symmetry>::operator [] (int i)
   { return elements_[i]; }
 
   /*
   * Element access operator (by reference).
   */
   template <class Symmetry> 
   inline 
   const Symmetry& SymmetryGroup<Symmetry>::operator [] (int i) const
   { return elements_[i]; }

   #ifndef PRDC_SYMMETRY_GROUP_TPP
   // Explicit instantiation declarations
   template <int D> class SpaceSymmetry;
   extern template class SymmetryGroup< SpaceSymmetry<1> >;
   extern template class SymmetryGroup< SpaceSymmetry<2> >;
   extern template class SymmetryGroup< SpaceSymmetry<3> >;
   #endif

}
}
#endif
