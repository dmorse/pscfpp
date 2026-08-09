#ifndef RP_FOURTH_ORDER_PARAMETER_BASE_TPP
#define RP_FOURTH_ORDER_PARAMETER_BASE_TPP

#include "FourthOrderParameterBase.h"

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/FFT.h>
#include <prdc/crystal/shiftToMinimum.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/Array.h>
#include <util/containers/DArray.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class T>
   FourthOrderParameterBase<D,T>::FourthOrderParameterBase(
                                       Simulator<D,T>& simulator,
                                       System<D,T>& system)
    : AverageAnalyzer<D,T>(simulator, system),
      kSize_(1),
      isInitialized_(false)
   {  ParamComposite::setClassName("FourthOrderParameter"); }

   /*
   * Setup before the main loop.
   */
   template <int D, class T>
   void FourthOrderParameterBase<D,T>::setup()
   {
      // Precondition: The system must have exactly two monomer types
      UTIL_CHECK(system().mixture().nMonomer() == 2);

      AverageAnalyzer<D,T>::setup();

      // Compute k-space mesh kMeshDimensions_ and kSize_
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      FFTT::computeKMesh(dimensions, kMeshDimensions_, kSize_);

      // Allocate variables
      if (!isInitialized_){
         wK_.allocate(dimensions);
         prefactor_.allocate(kMeshDimensions_);
         psi_.allocate(kMeshDimensions_);
         VecOp::eqS(prefactor_, 0.0);
      }

      computePrefactor();
      isInitialized_ = true;
   }

   /*
   * Compute and return the order parameter.
   */
   template <int D, class T>
   double FourthOrderParameterBase<D,T>::compute()
   {
      UTIL_CHECK(isInitialized_);
      UTIL_CHECK(wK_.capacity() == kSize_);
      UTIL_CHECK(prefactor_.capacity() == kSize_);
      UTIL_CHECK(psi_.capacity() == kSize_);
      UTIL_CHECK(system().w().hasData());

      if (!simulator().hasWc()){
         simulator().computeWc();
      }

      // Fourier transform W_(r) to obtain wK_ = W_(k)
      system().domain().fft().forwardTransform(simulator().wc(0), wK_);

      // Evaluate fourth powers, scaled by prefactors
      VecOp::sqSqAbsV(psi_, wK_);
      VecOp::mulEqV(psi_, prefactor_);

      // Summation
      double orderParameter = Reduce::sum(psi_, 1, kSize_);
      orderParameter = std::pow(orderParameter, 0.25);

      return orderParameter;
   }

   /*
   * Compute prefactors for all wavevectors.
   */
   template <int D, class T>
   void 
   FourthOrderParameterBase<D,T>::computePrefactor(Array<double>& prefactor)
   {
      IntVec<D> G;
      IntVec<D> Gmin;
      IntVec<D> nGmin;
      DArray< IntVec<D> > GminList;
      GminList.allocate(kSize_);
      MeshIterator<D> itr(kMeshDimensions_);
      MeshIterator<D> searchItr(kMeshDimensions_);

      // Calculate GminList
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      for (itr.begin(); !itr.atEnd(); ++itr){
         G = itr.position();
         Gmin = shiftToMinimum(G, meshDimensions, unitCell);
         GminList[itr.rank()] = Gmin;
      }

      // Compute prefactor for each wavevector
      for (itr.begin(); !itr.atEnd(); ++itr){
         bool inverseFound = false;

         // If prefactor of current wavevector has not been assigned
         if (prefactor[itr.rank()] == 0){
            Gmin = GminList[itr.rank()];

            // Compute inverse of wavevector
            nGmin.negate(Gmin);

            // Search for inverse of wavevector
            searchItr = itr;
            for (; !searchItr.atEnd(); ++searchItr){
               if (nGmin == GminList[searchItr.rank()]){
                  prefactor[itr.rank()] = 1.0/2.0;
                  prefactor[searchItr.rank()] = 1.0/2.0;
                  inverseFound = true;
               }
            }

            if (inverseFound == false){
               prefactor[itr.rank()]  = 1.0;
            }

         }

      }
   }

}
}
#endif
