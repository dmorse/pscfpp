#ifndef UTIL_SETABLE_H
#define UTIL_SETABLE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * Template for a value that can be set or declared null (i.e., unknown).
   *
   * Type T must be copy-constructable and have an assignment (=) operator.
   *
   * Convention for MPI programs: In parallel MPI programs in which a value 
   * for a variable is calculated by a reduce operation and is set only on 
   * a master processor, a default value should be set on all other processors 
   * whenever the true value is set on the master. This indicates on all 
   * processors that the value is known, though it may only be available
   * on the master processor. Similarly, when a value is unset, the unset()
   * function should be called on all processors. This convention allows
   * the isSet() function to be used on all processors to query whether the
   * value is known, which may be then be used to decide when to initiate
   * a recomputation that may require computation on all processors. This
   * convention is imposed by the isValid() function, which requires that
   * isSet have the same value on all processors within a communicator 
   * (i.e., all true or all false).
   *
   * \ingroup Misc_Module
   */
   template <class T>
   class Setable
   {

   public:

      /**
      * Default constructor.
      */
      Setable()
       : value_(),
         isSet_(false)
      {}

      /**
      * Copy constructor.
      *
      * \param other Setable object being copied.
      */
      Setable(const Setable<T>& other)
       : value_(other.value_),
         isSet_(other.isSet_)
      {}

      /**
      * Construct from T value (explicit).
      *
      * \param value  value of wrapped object
      */
      explicit Setable(const T& value)
       : value_(value),
         isSet_(true)
      {}

      /**
      * Assignment from another Setable<T> object.
      *
      * \param other  object on RHS of assignment
      */
      Setable<T>& operator = (const Setable<T>& other)
      {
         if (this != &other) {
            isSet_ = other.isSet_;
            if (other.isSet_) {
               value_ = other.value_;
            }
         }
         return *this;
      }

      /**
      * Assignment from T value.
      *
      * Equivalent to set(value). Sets the value and marks it as set.
      *
      * \param value T value on RHS of assignment
      * \return this object
      */
      Setable<T>& operator = (const T& value)
      {
         value_ = value;
         isSet_ = true;
         return *this;
      }

      /**
      * Set the value and mark as set.
      *
      * \param value  value to be assigned.
      */
      void set(const T& value)
      {
         value_ = value;
         isSet_ = true;
      }

      /**
      * Unset the value (mark as unknown).
      */
      void unset()
      {  isSet_ = false; }

      /**
      * Is this object set (is the value known)?
      *
      * \return true if set (known), false if null (unknown).
      */
      bool isSet() const
      {  return isSet_; }

      /**
      * Return value (if set).
      *
      * Throws an Exception if value is not set.
      */
      const T& value() const
      {
         if (!isSet_) {
            UTIL_THROW("Attempt to return unknown value.");
         }
         return value_;
      }

      #ifdef UTIL_MPI
      /**
      * Test consistency of states on different processors.
      *
      * If valid, return true, else throws an Exception.
      * The state is valid if the value of isSet is the
      * same on all processors.
      */
      bool isValid(MPI::Intracomm& communicator) const;
      #endif

   private:

      /// Value of associated variable.
      T value_;

      /// True if value is known (set), false if unknown (unset).
      bool isSet_;

   };

   #ifdef UTIL_MPI
   template <typename T>
   bool Setable<T>::isValid(MPI::Intracomm& communicator) const
   {
      int isSet = (int)isSet_;
      int total = 0;
      communicator.Allreduce(&isSet, &total, 1, MPI::INT, MPI::SUM);
      int nproc = communicator.Get_size();
      if (isSet_ && total != nproc) {
         UTIL_THROW("Inconsistent settings");
      }
      if ((!isSet_) && total != 0) {
         UTIL_THROW("Inconsistent settings");
      }
      return true;
   }
   #endif

}
#endif
