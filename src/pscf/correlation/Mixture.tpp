/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include "Polymer.h"

#include <pscf/chem/Composition.h>
#include <pscf/chem/SolventSpecies.h>

#include <util/global.h>

#include "Polymer.tpp"

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename WT>
   Mixture<WT>::Mixture()
    : mixturePtr_(nullptr)
   {}

   /*
   * Constructor.
   */
   template <typename WT>
   Mixture<WT>::Mixture(Composition<WT> const & mixture)
    : mixturePtr_(&mixture)
   {}

   /*
   * Destructor.
   */
   template <typename WT>
   Mixture<WT>::~Mixture()
   {}

   /*
   * Create an association with a Mixture.
   */
   template <typename WT>
   void Mixture<WT>::associate(Composition<WT> const & mixture)
   {
      UTIL_CHECK(!mixturePtr_);  
      mixturePtr_ = &mixture; 
   }

   /*
   * Allocate memory.
   */
   template <typename WT>
   void Mixture<WT>::allocate()
   {
      // Constant and preconditions
      UTIL_CHECK(mixturePtr_);  
      const int nMonomer = mixture().nMonomer();
      const int nPolymer = mixture().nPolymer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nPolymer > 0);
      UTIL_CHECK(!polymers_.isAllocated());

      // Allocation
      polymers_.allocate(nPolymer);
      for (int i = 0; i < nPolymer; ++i) {
         polymers_[i].associate(mixture().polymerSpecies(i));
         polymers_[i].allocate(nMonomer);
      }
   }

   /*
   * Compute mutable state variables by calling setup for all polymers.
   */
   template <typename WT>
   void Mixture<WT>::setup()
   {
      UTIL_CHECK(mixturePtr_);  
      UTIL_CHECK(isAllocated());  
      const int nMonomer = mixture().nMonomer();
      const int nPolymer = mixture().nPolymer();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(nPolymer > 0);
      UTIL_CHECK(polymers_.capacity() == nPolymer);

      // Construct array of segment lengths, indexed by monomer type
      DArray<double> kuhn;
      kuhn.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         kuhn[i] = mixture().monomer(i).kuhn();
      }

      // Setup all Correlation::Polymer objects
      for (int i = 0; i < nPolymer; ++i) {
         polymers_[i].setup(kuhn);
      }

   }

   /*
   * Compute an array of Omega(k) values for one monomer type pair.
   */
   template <typename WT>
   void Mixture<WT>::computeOmega(int ma, int mb, 
                              Array<double> const & kSq, 
                              Array<double> & correlations) const
   {
      // Constants
      const double vMonomer = mixture().vMonomer();
      const int nPolymer = mixture().nPolymer();
      const int nSolvent = mixture().nSolvent();
      const int nk = kSq.capacity();

      // Preconditions
      UTIL_CHECK(nPolymer > 0);
      UTIL_CHECK(polymers_.capacity() == nPolymer);
      UTIL_CHECK(nSolvent >= 0);
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(correlations.capacity() == nk);

      // Initialize correlations array to zero
      for (int j = 0; j < nk; ++j){
         correlations[j] = 0.0;
      }

      // Loop over polymer species
      double phi, length, c;
      int ip;      // polymer species index
      int nA, nB;  // number of blocks of types ma and mb 
      int ja, jb;  // block indices
      int ia, ib;  // loop indices
      for (ip = 0; ip < nPolymer; ip++){

         // Polymer properties
         Polymer<WT> const & polymer = polymers_[ip];
         phi = polymer.phi();
         length = polymer.totalLength();
         c = phi/(length*vMonomer);

         if (ma == mb) {

            // Case: Equal monomer type indices
            nA = polymer.blockIds(ma).size(); // # of blocks of type ma
            if (nA > 0) {
               // Intrablock contribution
               for (ia = 0; ia < nA; ++ia) {
                   ja = polymer.blockIds(ma)[ia]; // block index
                   polymer.computeOmega(ja, ja, c, kSq, correlations);
               }
               if (nA > 1) {
                  // If the polymer contains multiple blocks of type ma,
                  // compute inter-block contributions
                  for (ia = 0; ia < nA; ++ia) {
                     ja = polymer.blockIds(ma)[ia]; // block index
                     for (ib = 0; ib < ia; ++ib) {
                        jb = polymer.blockIds(ma)[ib]; // block index
                        UTIL_CHECK(jb != ja);
                        polymer.computeOmega(ja, jb, 2.0*c, kSq, 
                                             correlations);
                     }
                  }
               }
            }

         } else {

            // Case: Unequal monomer type indices
            nA = polymer.blockIds(ma).size(); // # of blocks of type ma
            nB = polymer.blockIds(mb).size(); // # of blocks of type mb
            if (nA > 0 && nB > 0) {
               for (ia = 0; ia < nA; ++ia) {
                  ja = polymer.blockIds(ma)[ia]; // block index for ma
                  for (ib = 0; ib < nB; ++ib) {
                     jb = polymer.blockIds(mb)[ib]; // block index for  mb
                     UTIL_CHECK(jb != ja);
                     polymer.computeOmega(ja, jb, c, kSq, correlations);
                  }
               }
            }

         }

      } // end loop over polymer species

      // Solvent contributions (if any)
      if (nSolvent > 0 and ma == mb) {
         double size;
         for (int i = 0; i < nSolvent; i++) {
            SolventSpecies<WT> const & solvent 
                                         = mixture().solventSpecies(i);
            if (solvent.monomerId() == ma) {
               phi = solvent.phi();
               size = solvent.size();
               c = phi*size/vMonomer;
               for (int i = 0; i < nk; ++i) {
                  correlations[i] += c;
               }
            }
         }
      }

   }

   /*
   * Compute k-space array of total intramolecular correlation function.
   */
   template <typename WT>
   void Mixture<WT>::computeOmegaTotal(Array<double> const & kSq, 
                                       Array<double> & correlations) const
   {
      // Constants
      const double vMonomer = mixture().vMonomer();
      const int nPolymer = mixture().nPolymer();
      const int nSolvent = mixture().nSolvent();
      const int nk = kSq.capacity();

      // Preconditions
      UTIL_CHECK(nPolymer > 0);
      UTIL_CHECK(nSolvent >= 0);
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(correlations.capacity() == nk);

      // Initialize correlations to zero
      for (int j = 0; j < nk; ++j){
         correlations[j] = 0.0;
      }

      // Loop over polymer species
      double phi, totalLength, c;
      for (int ip = 0; ip < nPolymer; ip++){
         Correlation::Polymer<WT> const & polymer = polymers_[ip];
         phi = polymer.phi();
         totalLength = polymer.totalLength();
         c = phi/(totalLength*vMonomer);
         polymer.computeOmegaTotal(c, kSq, correlations);
      }

      // Loop over solvent species (if any)
      if (nSolvent > 0) {
         double size;
         for (int i = 0; i < nSolvent; i++){
            SolventSpecies<WT> const & solvent 
                                     = mixture().solventSpecies(i);
            phi = solvent.phi();
            size = solvent.size();
            c = phi*size/vMonomer;
            for (int j = 0; j < nk; ++j) {
               correlations[j] += c;
            }
         }
      }

   }

}
}
