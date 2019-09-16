#ifndef PSCF_SYMMETRY_GROUP_H
#define PSCF_SYMMETRY_GROUP_H

/*
* Pscfatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

//#include <iostream>
#include <vector>

namespace Pscf 
{

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
   * \ingroup Crystal_Module
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
      * Assignment operator.
      */
      SymmetryGroup<Symmetry>& 
      operator = (const SymmetryGroup<Symmetry>& other);
 
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


   // Member function definitions (non-inline)

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

   #ifndef PSCF_SYMMETRY_GROUP_TPP
   template <int D> class SpaceSymmetry;
   extern template class SymmetryGroup< SpaceSymmetry<1> >;
   extern template class SymmetryGroup< SpaceSymmetry<2> >;
   extern template class SymmetryGroup< SpaceSymmetry<3> >;
   #endif 

}
// #include "SymmetryGroup.tpp"
#endif
