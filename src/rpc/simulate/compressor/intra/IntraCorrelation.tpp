#ifndef RPC_INTRACORRELATION_TPP
#define RPC_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"
#include <rpc/System.h>
#include <util/global.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/PolymerType.h>
#include <prdc/crystal/shiftToMinimum.h>

namespace Pscf {
namespace Rpc{

   using namespace Util;

   // Constructor
   template <int D>
   IntraCorrelation<D>::IntraCorrelation(System<D>& system)
    : systemPtr_(&system),
      kSize_(1)
   {  setClassName("IntraCorrelation"); }

   // Destructor
   template <int D>
   IntraCorrelation<D>::~IntraCorrelation()
   {}

   template<int D>
   double IntraCorrelation<D>::computeDebye(double x)
   {
      if (x == 0){
         return 1.0;
      } else {
         return 2.0 * (std::exp(-x) - 1.0 + x) / (x * x);
      }
   }

   template<int D>
   double IntraCorrelation<D>::computeIntraCorrelation(double qSquare)
   {
      const int np = system().mixture().nPolymer();
      const double vMonomer = system().mixture().vMonomer();
      
      // Overall intramolecular correlation
      double omega = 0;
      int monomerId; int nBlock;
      double phi; double kuhn; double length; 
      double totalN; double avgKuhn; double g; double rg2;
      Polymer<D> const * polymerPtr;
      PolymerType::Enum type;
      for (int i = 0; i < np; i++){
         polymerPtr = &system().mixture().polymer(i);
         
         // Ensure this function only works for linear polymers
         type = polymerPtr->type();
         UTIL_CHECK(type == PolymerType::Linear);
         
         phi = polymerPtr->phi();
         nBlock = polymerPtr->nBlock();
         totalN = 0;
         avgKuhn = 0;
         for (int j = 0; j < nBlock; j++) {
            monomerId = polymerPtr-> block(j).monomerId();
            kuhn = system().mixture().monomer(monomerId).kuhn();
            length = polymerPtr-> block(j).length();
            totalN += length;
            avgKuhn += kuhn/nBlock;
         }
         rg2 = totalN* avgKuhn* avgKuhn /6.0;
         g = computeDebye(qSquare*rg2);
         omega += phi*totalN*g/ vMonomer;
      }
      
      return omega;
   }

   template<int D>
   RField<D> IntraCorrelation<D>::computeIntraCorrelations()
   {
      
      IntVec<D> const & dimensions = system().mesh().dimensions();
      
      // Compute Fourier space kMeshDimensions_
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = dimensions[i];
            kSize_ *= dimensions[i];
         } else {
            kMeshDimensions_[i] = dimensions[i]/2 + 1;
            kSize_ *= (dimensions[i]/2 + 1);
         }
      }
      
      // Define IntraCorrelations.
      RField<D> intraCorrelations;
      intraCorrelations.allocate(kMeshDimensions_);
      
      // Define iterator
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, system().mesh().dimensions(), 
                                  system().unitCell());
         Gsq = system().unitCell().ksq(Gmin);
         intraCorrelations[iter.rank()] = computeIntraCorrelation(Gsq);
      }
      
      return intraCorrelations;
   }

}
}
#endif
