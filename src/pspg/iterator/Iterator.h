#ifndef PSPG_ITERATOR_H
#define PSPG_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <prdc/cuda/Field.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   typedef Field<cudaReal> FieldCUDA;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspg_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite
   {

   public:

      #if 0
      /**
      * Default constructor.
      */
      Iterator();
      #endif

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
      bool isSymmetric() const
      {  return (isSymmetric_); }

      /**
      * Is the unit cell flexible (true) or rigid (false).
      */
      bool isFlexible() const;

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
      * Return reference to parent system.
      */
      System<D>& system();

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;
      
   };

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

   template<int D>
   inline bool Iterator<D>::isFlexible() const
   {  return isFlexible_; }

   template<int D>
   inline System<D>& Iterator<D>::system() 
   {  return *sysPtr_; }

} // namespace Pspg
} // namespace Pscf
#endif
