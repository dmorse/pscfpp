#ifndef RP_ITERATOR_H
#define RP_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/FSArray.h>      // member
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Base class template for iterative solvers for SCF equations.
   *
   * Template parameters:
   *
   *    - D  dimension of space
   *    - T  Types class (Rpc::Types<D> or Rpg::Types<D>)
   *
   * \ingroup Rp_Scft_Module
   */
   template <int D, class T>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Iterate to solution.
      *
      * \param isContinuation true iff a continuation within a sweep
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve(bool isContinuation) = 0;

      /**
      * Log output timing results
      */
      virtual void outputTimers(std::ostream& out) const = 0;

      /**
      * Clear timers
      */
      virtual void clearTimers() = 0;

      /**
      * Does this iterator use a symmetry-adapted Fourier basis?
      */
      bool isSymmetric() const
      {  return (isSymmetric_); }

      /**
      * Return true iff unit cell has any flexible lattice parameters.
      */
      bool isFlexible() const
      {  return (isFlexible_); }

      /**
      * Get the array indicating which lattice parameters are flexible.
      *
      * This array should be nParameters long, where the i-th entry is a
      * boolean indicating whether parameter i is flexible.
      */
      FSArray<bool,6> flexibleParams() const
      {  return flexibleParams_; }

      /**
      * Get the number of flexible lattice parameters.
      */
      int nFlexibleParams() const;

      /**
      * Set the array indicating which lattice parameters are flexible.
      *
      * \param flexParams input boolean array
      */
      void setFlexibleParams(FSArray<bool,6> const & flexParams);

      /**
      * Return the stress used by this Iterator, for one lattice parameter.
      *
      * Will throw an error if paramId corresponds to a lattice parameter
      * that is not flexible (according to the flexibleParams array).
      *
      * \param paramId  index of lattice parameter
      */
      virtual double stress(int paramId) const;

   protected:

      /**
      * Does this iterator use a symmetry-adapted basis?
      */
      bool isSymmetric_;

      /**
      * Are any lattice parameters flexible during iteration?
      */
      bool isFlexible_;

      /**
      * Array of indices of the lattice parameters that are flexible.
      */
      FSArray<bool,6> flexibleParams_;

      // Protected member functions

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      Iterator(typename T::System& system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Set parent system.
      *
      * \param system  parent System object
      */
      void setSystem(typename T::System& system);

      /**
      * Get parent system by const reference.
      */
      typename T::System const & system() const;

      /**
      * Get parent system by non-const reference.
      */
      typename T::System& system();

   private:

      /// Pointer to the associated system object.
      typename T::System* sysPtr_;

   };

   // Inline member functions

   /**
   * Get parent system by const reference.
   */
   template <int D, class T> inline
   typename T::System const & Iterator<D, T>::system() const
   {
      UTIL_CHECK(sysPtr_);
      return *sysPtr_;
   }

   /**
   * Get parent system by non-const reference.
   */
   template <int D, class T> inline
   typename T::System& Iterator<D, T>::system()
   {
      UTIL_CHECK(sysPtr_);
      return *sysPtr_;
   }

} // namespace Rp
} // namespace Pscf
#endif
