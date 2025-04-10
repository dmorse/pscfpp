#ifndef RPG_INTRACORRELATION_TPP
#define RPG_INTRACORRELATION_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"
#include <rpg/System.h>
#include <util/global.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/chem/Debye.h>
#include <pscf/chem/EdgeIterator.h>
#include <pscf/iterator/NanException.h>
#include <pscf/chem/PolymerType.h>
#include <pscf/cuda/HostDArray.h>
#include <prdc/crystal/shiftToMinimum.h>

namespace Pscf {
namespace Rpg{

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
   void IntraCorrelation<D>::computeIntraCorrelations(RField<D>& intraCorrelations)
   {
      const double vMonomer = system().mixture().vMonomer();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      
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
      
      // CudaReal array
      HostDArray<cudaReal> hostIntra;
      hostIntra.allocate(kSize_);
      
      // Initialize intraCorrelations
      for (int i = 0; i < hostIntra.capacity(); ++i){
         hostIntra[i] = 0.0;
      }

      // Overall intramolecular correlation
      int nPolymer = system().mixture().nPolymer();
      int nBlock; int monomerId;
      double phi; double kuhn; double length; double totalN; 
      double prefactor; double ksq; 
      int monomerIdA; int monomerIdB;
      double eA; double eB; double lengthA; double lengthB;
      double kuhnA; double kuhnB;

      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      DArray<double> Gsq;
      Gsq.allocate(kSize_);
      UnitCell<D> const & unitCell = system().domain().unitCell();
      
      // Compute qSq
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, dimensions, unitCell);
         Gsq[iter.rank()] = unitCell.ksq(Gmin);
      }
      
      // Loop over polymer species  
      for (int i = 0; i < nPolymer; i++){
         Mixture<D> const & mixture = system().mixture();
         Polymer<D> const & polymer = mixture.polymer(i);
         
         phi = polymer.phi();
         totalN = 0;
         nBlock = polymer.nBlock();
         
         // Compute totalN 
         for (int j = 0; j < nBlock; j++) {
            length = polymer.block(j).length();
            totalN += length;
         }
         
         prefactor = phi/(totalN*vMonomer);
            
         // Compute diagonal (j = k) contributions
         for (int j = 0; j < nBlock; j++) {
            monomerId = polymer.block(j).monomerId();
            kuhn = mixture.monomer(monomerId).kuhn();
            length = polymer.block(j).length();
            
            // Loop over ksq
            for (iter.begin(); !iter.atEnd(); ++iter) {
               ksq = Gsq[iter.rank()];
               hostIntra[iter.rank()] += prefactor * Debye::d(ksq, length, kuhn);
            }
            
         }
         
         // Compute diagonal contributions
         if (nBlock > 1) {
            prefactor = 2.0*prefactor;
            EdgeIterator EdgeItr(polymer);
            
            for (int ia = 1; ia < nBlock; ++ia) {
               Block<D> const & blockA = polymer.block(ia);
               lengthA = blockA.length();
               monomerIdA = blockA.monomerId();
               kuhnA = mixture.monomer(monomerIdA).kuhn();
               for (int ib = 0; ib < ia; ++ib)  {
                  Block<D> const & blockB = polymer.block(ib);
                  lengthB = blockB.length();
                  monomerIdB = blockB.monomerId();
                  kuhnB = mixture.monomer(monomerIdB).kuhn();
                  
                  // Find path from edge ia to ib
                  EdgeItr.begin(ia, ib);
                  
                  // Mean-square end-to-end length of segment between two blocks
                  double d = 0;

                  while (EdgeItr.notEnd()) {
                     if (EdgeItr.currentEdgeId() != ia && EdgeItr.currentEdgeId() != ib){
                        Block<D> const & blockK = polymer.block(EdgeItr.currentEdgeId());
                        double lengthK = blockK.length();
                        double monomerIdK = blockK.monomerId();
                        double kuhnK = mixture.monomer(monomerIdK).kuhn();
                        d += lengthK * kuhnK * kuhnK;
                     }
                     ++EdgeItr;
                  }
                  
                  double x;
                  
                  // Loop over ksq
                  for (iter.begin(); !iter.atEnd(); ++iter) {
                     ksq = Gsq[iter.rank()];
                     x = d * ksq/6.0;
                     eA = Debye::e(ksq, lengthA, kuhnA);
                     eB = Debye::e(ksq, lengthB, kuhnB);
                     hostIntra[iter.rank()] += std::exp(-x) * prefactor * eA * eB;
                  }
               }
            }
         }
      }
      
      // Copy to device (gpu) memory
      intraCorrelations = hostIntra;
      
   }

}
}
#endif
