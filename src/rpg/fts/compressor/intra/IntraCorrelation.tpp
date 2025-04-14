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
   void 
   IntraCorrelation<D>::computeIntraCorrelations(RField<D>& correlations)
   {
      // Local copies of system properties
      Mixture<D> const & mixture = system().mixture();
      UnitCell<D> const & unitCell = system().domain().unitCell();
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      const int nPolymer = mixture.nPolymer();
      const int nSolvent = mixture.nSolvent();
      const double vMonomer = mixture.vMonomer();

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
      UTIL_CHECK(correlations.capacity() == kSize_);

      // Allocate Gsq (k-space array of square wavenumber values)
      DArray<double> Gsq;
      Gsq.allocate(kSize_);

      // Compute Gsq
      IntVec<D> G, Gmin;
      MeshIterator<D> iter;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         G = iter.position();
         Gmin = shiftToMinimum(G, dimensions, unitCell);
         Gsq[iter.rank()] = unitCell.ksq(Gmin);
      }

      double phi, cPolymer, polymerLength;
      double length, lengthA, lengthB, ksq;
      double kuhn, kuhnA, kuhnB, eA, eB;
      int monomerId, monomerIdA, monomerIdB, rank;

      // Allocate local array
      HostDArray<cudaReal> correlations_h;
      correlations_h.allocate(kSize_);
      
      // Initialize correlations_h to zero
      for (int i = 0; i < correlations_h.capacity(); ++i){
         correlations_h[i] = 0.0;
      }

      // Loop over polymer species
      for (int i = 0; i < nPolymer; i++){

         // Local copies of polymer properties
         Polymer<D> const & polymer = mixture.polymer(i);
         const int nBlock = polymer.nBlock();
         phi = polymer.phi();

         // Compute polymerLength (sum of lengths of all blocks)
         polymerLength = 0.0;
         for (int j = 0; j < nBlock; j++) {
            length = polymer.block(j).length();
            polymerLength += length;
         }

         // Compute cPolymer (polymer number concentration)
         cPolymer = phi/(polymerLength*vMonomer);

         // Compute diagonal (j = k) contributions
         for (int j = 0; j < nBlock; j++) {

            monomerId = polymer.block(j).monomerId();
            kuhn = mixture.monomer(monomerId).kuhn();
            length = polymer.block(j).length();

            // Loop over ksq to increment correlations_h
            for (iter.begin(); !iter.atEnd(); ++iter) {
               rank = iter.rank();
               ksq = Gsq[rank];
               correlations_h[rank] += cPolymer * Debye::d(ksq, length, kuhn);
            }

         }

         // Compute off-diagonal contributions
         if (nBlock > 1) {
            EdgeIterator EdgeItr(polymer);

            // Outer loop over blocks
            for (int ia = 1; ia < nBlock; ++ia) {

               // Block A properties
               Block<D> const & blockA = polymer.block(ia);
               lengthA = blockA.length();
               monomerIdA = blockA.monomerId();
               kuhnA = mixture.monomer(monomerIdA).kuhn();

               // Inner loop over blocks
               for (int ib = 0; ib < ia; ++ib)  {

                  // Block B properties
                  Block<D> const & blockB = polymer.block(ib);
                  lengthB = blockB.length();
                  monomerIdB = blockB.monomerId();
                  kuhnB = mixture.monomer(monomerIdB).kuhn();

                  // Mean-square end-to-end length of segment between blocks
                  double d = 0.0;
                  int edgeId;
                  EdgeItr.begin(ia, ib);
                  while (EdgeItr.notEnd()) {
                     edgeId = EdgeItr.currentEdgeId();
                     if (edgeId != ia && edgeId != ib){
                        Block<D> const & blockK = polymer.block(edgeId);
                        double lengthK = blockK.length();
                        double monomerIdK = blockK.monomerId();
                        double kuhnK = mixture.monomer(monomerIdK).kuhn();
                        d += lengthK * kuhnK * kuhnK;
                     }
                     ++EdgeItr;
                  }

                  // Loop over ksq to increment correlations_h
                  double x;
                  double prefactor = 2.0*cPolymer;
                  for (iter.begin(); !iter.atEnd(); ++iter) {
                     rank = iter.rank();
                     ksq = Gsq[rank];
                     x = std::exp( -d * ksq / 6.0);
                     eA = Debye::e(ksq, lengthA, kuhnA);
                     eB = Debye::e(ksq, lengthB, kuhnB);
                     correlations_h[rank] += prefactor * x * eA * eB;
                  }

               } // loop: block ib
            } // loop: block ia
         } // if (nBlock > 1)

      } // loop over polymers


      // Loop over solvent species (if any)
      if (nSolvent > 0) {
         double size, dcorr;
         for (int i = 0; i < nSolvent; i++){
            Solvent<D> const & solvent = mixture.solvent(i);
            phi = solvent.phi();
            size = solvent.size();
            dcorr = phi*size/vMonomer;
            for (iter.begin(); !iter.atEnd(); ++iter) {
               correlations_h[iter.rank()] += dcorr;
            }
         }
      }

      // Copy local array correlations_h to device (gpu) memory
      correlations = correlations_h;

   }

}
}
#endif
