#ifndef RPG_ITERATOR_H
#define RPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class

#include <prdc/cuda/types.h>              // typedef
#include <pscf/cuda/DeviceArray.h>        // typedef
#include <util/containers/FSArray.h>      // member

namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   typedef DeviceArray<cudaReal> FieldCUDA;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      Iterator(System<D>& system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Iterate to solution.
      *
      * \param isContinuation  true iff continuation within a sweep
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
      bool isSymmetric() const;

      /**
      * Is the unit cell flexible (true) or rigid (false).
      */
      bool isFlexible() const;

      /**
      * Get the array indicating which lattice parameters are flexible.
      *
      * This array should be nParameters long, where the i-th entry is a
      * boolean indicating whether parameter i is flexible.
      */
      FSArray<bool, 6> flexibleParams() const;

      /**
      * Get the number of flexible lattice parameters.
      */
      int nFlexibleParams() const;

      /**
      * Set the array indicating which lattice parameters are flexible.
      *
      * \param flexParams input boolean array
      */
      void setFlexibleParams(FSArray<bool, 6> const & flexParams);

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
      * Is the unit cell flexible during iteration?
      */
      bool isFlexible_;

      /**
      * Array of indices of the lattice parameters that are flexible.
      */
      FSArray<bool, 6> flexibleParams_;

      /**
      * Return reference to parent system.
      */
      System<D>& system();

      /**
      * Return const reference to parent system.
      */
      System<D> const & system() const;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Protected inline member functions

   // Return reference to parent system.
   template<int D> inline 
   System<D>& Iterator<D>::system()
   {  return *sysPtr_; }

   // Return const reference to parent system.
   template<int D> inline 
   System<D> const & Iterator<D>::system() const
   {  return *sysPtr_; }

   // Explicit instantiation declarations
   extern template class Iterator<1>;
   extern template class Iterator<2>;
   extern template class Iterator<3>;

} // namespace Rpg
} // namespace Pscf
#endif
