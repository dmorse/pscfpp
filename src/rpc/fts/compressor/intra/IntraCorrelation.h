#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>      // memmber variable type

// Forward references
namespace Pscf {
   namespace Prdc {
      namespace Cpu {
         template <int D> class RField;
      }
   }
   namespace Rpc {
      template <int D> class System;
   }
}


namespace Pscf {
namespace Rpc {

   using namespace Pscf::Prdc::Cpu;

   /**
   * Intramolecular correlation analysis for LR compressors.
   *
   * \ingroup Rpc_Compressor_Intra_Module
   */
   template <int D>
   class IntraCorrelation 
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
      * Compute and modify intramolecular correlations.
      *
      * \param correlations  k-space grid of omega values
      */
      void computeIntraCorrelations(RField<D>& correlations);

   protected:

      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
   private:
      
      /// Pointer to the associated system object.
      System<D>* systemPtr_;
      
      /// Dimensions of Fourier grid for DFT of a real function
      IntVec<D> kMeshDimensions_;
      
      /// Number of elements in the Fourier grid
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
