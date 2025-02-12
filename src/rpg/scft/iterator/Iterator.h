#ifndef RPG_ITERATOR_H
#define RPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpg/scft/sweep/Sweep.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/sweep/ParameterModifier.h> // base class
#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>                  

namespace Pscf {
namespace Rpg
{

   template <int D>
   class System;

   using namespace Util;

   typedef Prdc::Cuda::DeviceArray<Prdc::Cuda::cudaReal> FieldCUDA;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Rpg_Scft_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite, public ParameterModifier
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
      virtual void outputTimers(std::ostream& out) = 0;
      
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
      FSArray<bool,6> flexibleParams() const;

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
      FSArray<bool,6> flexibleParams_;

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

   // Default constructor
   template<int D>
   inline Iterator<D>::Iterator()
   : isSymmetric_(false),
     isFlexible_(false),
     sysPtr_(nullptr)
   {  setClassName("Iterator"); }

   // Constructor
   template<int D>
   inline Iterator<D>::Iterator(System<D>& system)
   : isSymmetric_(false),
     isFlexible_(false),
     sysPtr_(&system)
   {  setClassName("Iterator"); }

   // Destructor
   template<int D>
   inline Iterator<D>::~Iterator()
   {}

   // Does this iterator use a symmetry-adapted Fourier basis?
   template<int D>
   inline bool Iterator<D>::isSymmetric() const
   {  return isSymmetric_; }


   // Is the unit cell flexible (true) or rigid (false) ?
   template<int D>
   inline bool Iterator<D>::isFlexible() const
   {  return isFlexible_; }

   // Get the array indicating which lattice parameters are flexible.
   template<int D>
   inline FSArray<bool,6> Iterator<D>::flexibleParams() const
   {  return flexibleParams_; }

   // Get the number of flexible lattice parameters
   template <int D>
   int Iterator<D>::nFlexibleParams() const
   {
      UTIL_CHECK(flexibleParams_.size() == 
                                 system().unitCell().nParameter());
      int nFlexParams = 0;
      for (int i = 0; i < flexibleParams_.size(); i++) {
         if (flexibleParams_[i]) nFlexParams++;
      }
      return nFlexParams;
   }

   // Set the array indicating which lattice parameters are flexible.
   template <int D>
   void Iterator<D>::setFlexibleParams(FSArray<bool,6> const & flexParams)
   {  
      flexibleParams_ = flexParams; 
      if (nFlexibleParams() == 0) {
         isFlexible_ = false;
      } else {
         isFlexible_ = true;
      }
   }

   // Return reference to parent system.
   template<int D>
   inline System<D>& Iterator<D>::system() 
   {  return *sysPtr_; }

   // Return const reference to parent system.
   template<int D>
   inline System<D> const & Iterator<D>::system() const
   {  return *sysPtr_; }

} // namespace Rpg
} // namespace Pscf
#endif
