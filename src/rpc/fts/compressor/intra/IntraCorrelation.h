#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <prdc/cpu/RField.h>              // member
#include <prdc/cpu/RFieldDft.h>           // member
#include <util/containers/DArray.h>       // member

namespace Pscf {
namespace Rpc
{

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Base class for iterators that impose incompressibility.
   *
   * \ingroup Rpc_Compressor_Intra_Module
   */
   template <int D>
   class IntraCorrelation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object
      */
      IntraCorrelation(System<D>& system);

      /**
      * Destructor.
      */
      ~IntraCorrelation();

      /**
      * Compute and modify intramolecular correlations 
      */
      void computeIntraCorrelations(RField<D>& intraCorrelations);

   protected:

      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
   private:
      
      /// Pointer to the associated system object.
      System<D>* systemPtr_;
      
      /// Dimensions of fourier space 
      IntVec<D> kMeshDimensions_;
      
      /// Size of fourier space
      int kSize_;
      
   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& IntraCorrelation<D>::system()
   {  return *systemPtr_; }

   
   #ifndef RPC_INTRACORRELATION_TPP
   // Suppress implicit instantiation
   extern template class IntraCorrelation<1>;
   extern template class IntraCorrelation<2>;
   extern template class IntraCorrelation<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif
