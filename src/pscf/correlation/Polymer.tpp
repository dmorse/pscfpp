#ifndef PSCF_CORRELATION_POLYMER_TPP
#define PSCF_CORRELATION_POLYMER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/Edge.h>
#include <pscf/chem/EdgeIterator.h>
#include <pscf/chem/PolymerModel.h>
#include <pscf/correlation/Debye.h>
#include <util/global.h>

#include <cmath>

namespace Pscf {
namespace Correlation {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename WT>
   Polymer<WT>::Polymer()
    : speciesPtr_(nullptr),
      phi_(-1.0),
      totalLength_(-1.0),
      nBlock_(0)
   {}

   /*
   * Constructor.
   */
   template <typename WT>
   Polymer<WT>::Polymer(PolymerSpecies<WT> const & polymer)
    : speciesPtr_(&polymer),
      phi_(-1.0),
      totalLength_(-1.0),
      nBlock_(0)
   {}

   /*
   * Destructor.
   */
   template <typename WT>
   Polymer<WT>::~Polymer()
   {}

   /*
   * Create an association with a PolymerSpecies descriptor.
   */
   template <typename WT>
   void Polymer<WT>::associate(PolymerSpecies<WT> const & polymer)
   {  speciesPtr_ = &polymer; }

   /*
   * Allocation and initialization that depend on immutable data.
   */
   template <typename WT>
   void Polymer<WT>::allocate(int nMonomer)
   {
      UTIL_CHECK(speciesPtr_);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
      nBlock_ = species().nBlock();
      UTIL_CHECK(nBlock_ > 0);

      kuhn_.allocate(nBlock_);
      length_.allocate(nBlock_);
      rSq_.allocate(nBlock_, nBlock_);
      blockIds_.allocate(nMonomer_);

      // Construct blocksIds_ container
      int monomerId;
      for (int i = 0; i < nBlock_; ++i) {
         monomerId = species().edge(i).monomerId();
         UTIL_CHECK(monomerId >= 0);
         UTIL_CHECK(monomerId < nMonomer_);
         blockIds_[monomerId].append(i);
      }
   }

   /*
   * Set member variables that depends on mutable system data.
   */
   template <typename WT>
   void Polymer<WT>::setup(Array<double> const & kuhn)
   {
      // Preconditions
      UTIL_CHECK(speciesPtr_);
      UTIL_CHECK(nBlock_ > 0);
      UTIL_CHECK(kuhn.capacity() == nMonomer_);

      /*
      * Note: Array parameter kuhn contains statistical segment lengths
      * indexed by monomer type id. Member variable kuhn_ contains segment
      * length values for specific blocks, indexed by block id.
      */

      // Set block properties (arrays length_ and kuhn_) and totalLength_
      int monomerId;
      totalLength_ = 0.0;
      for (int i = 0; i < nBlock_; ++i) {
         Edge const & edge = species().edge(i);
         if (PolymerModel::isThread()) {
            length_[i] = edge.length();
         } else {
            length_[i] = (double) edge.nBead();
         }
         totalLength_ += length_[i];
         monomerId = edge.monomerId();
         // std::cout << "\n i, momonomerId, length = " 
         //          << i << "  " <<  monomerId << " " << length_[i];
         UTIL_CHECK(monomerId >= 0);
         UTIL_CHECK(monomerId < nMonomer_);
         kuhn_[i] = kuhn[monomerId];
      }
      // std::cout << "\n totalLength = " << totalLength_;

      // Set species volume fraction phi_
      phi_ = species().phi();

      // Set all diagonal elements of rSq_ matrix to zero
      for (int ia = 0; ia < nBlock_; ++ia) {
         rSq_(ia, ia) = 0.0; 
      }

      // Compute off-diagonal elements of Rsq_
      if (nBlock_ > 1) {

         EdgeIterator<WT> edgeItr(species());
         double kuhnA, kuhnB, kuhnC, lengthC, rsq;
         int ia, ib, ic;

         // Outer loop over block A
         for (ia = 1; ia < nBlock_; ++ia) {
            kuhnA = kuhn_[ia];
   
            // Inner loop over block B (ib < ia)
            for (ib = 0; ib < ia; ++ib)  {
               kuhnB = kuhn_[ib];
   
               // Initialize rsq
               if (PolymerModel::isThread()) {
                  rsq = 0.0;
               } else {
                  rsq = 0.5*(kuhnA*kuhnA + kuhnB*kuhnB);
               }
   
               // Loop over intermediate block C, if any
               edgeItr.begin(ia, ib);
               while (edgeItr.notEnd()) {
                  ic = edgeItr.currentEdgeId();
                  if (ic != ia && ic != ib) {
                     lengthC = length_[ic];
                     kuhnC = kuhn_[ic];
                     rsq += lengthC * kuhnC * kuhnC;
                     // Includes two half-bonds at ends of block C
                  }
                  ++edgeItr;
               }
   
               // Assign matrix elements of rSq_
               rSq_(ia, ib) = rsq;
               rSq_(ib, ia) = rsq;
   
            } // end loop over block B
         } // end loop over block A
      } // end if nBlock_ > 1 

   }

   /*
   * Compute and return an intramolecular correlation functions.
   */
   template <typename WT>
   double Polymer<WT>::computeOmega(int ia, int ib, 
                                double prefactor, 
                                double kSq) const
   {
      // Preconditions
      UTIL_CHECK(speciesPtr_);         // check if associated
      UTIL_CHECK(nBlock_ > 0);         // check if allocated
      UTIL_CHECK(totalLength_ > 0.0);  // check if setup

      double correlation;
      double lengthA = length_[ia];
      double kuhnA = kuhn_[ia];
      if (ia == ib) {
         if (PolymerModel::isThread()) {
            correlation = prefactor * Correlation::dt(kSq, lengthA, kuhnA);
         } else {
            correlation = prefactor * Correlation::db(kSq, lengthA, kuhnA);
         }
      } else {
         double lengthB = length_[ib];
         double kuhnB = kuhn_[ib];
         double x = std::exp( -rSq_(ia, ib) * kSq / 6.0);
         double eA, eB;
         if (PolymerModel::isThread()) {
            eA = Correlation::et(kSq, lengthA, kuhnA);
            eB = Correlation::et(kSq, lengthB, kuhnB);
         } else {
            eA = Correlation::eb(kSq, lengthA, kuhnA);
            eB = Correlation::eb(kSq, lengthB, kuhnB);
         }
         correlation = prefactor * x * eA * eB;
      }
      return correlation;
   }

   /*
   * Compute array of correlation function values for a pair of blocks.
   *
   * Parameters:
   *   - ia and ib are block indices.
   *   - kSq and correlation must have the same nonzero capacity
   */
   template <typename WT>
   void Polymer<WT>::computeOmega(int ia, int ib, 
                              double prefactor, 
                              Array<double> const & kSq,
                              Array<double> & correlation) const
   {
      // Preconditions
      UTIL_CHECK(speciesPtr_);          // check if associated
      UTIL_CHECK(nBlock_ > 0);          // check if allocated
      UTIL_CHECK(totalLength_ > 0.0);   // check if setup
      const int nk = kSq.capacity();
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(correlation.capacity() == nk);

      const double lengthA = length_[ia];
      const double kuhnA = kuhn_[ia];

      if (ia == ib) {
         double d;
         if (PolymerModel::isThread()) {
            for (int j = 0; j < nk; ++j) {
               d = Correlation::dt(kSq[j], lengthA, kuhnA);
               correlation[j] += prefactor * d;
            }
         } else {
            for (int j = 0; j < nk; ++j) {
               d = Correlation::db(kSq[j], lengthA, kuhnA);
               correlation[j] += prefactor * d;
            }
         }
      } else {
         const double lengthB = length_[ib];
         const double kuhnB = kuhn_[ib];
         const double alpha = -1.0*rSq_(ia, ib)/6.0;
         double ks, x, eA, eB;
         if (PolymerModel::isThread()) {
            for (int j = 0; j < nk; ++j) {
               ks = kSq[j];
               x = std::exp(alpha*ks);
               eA = Correlation::et(ks, lengthA, kuhnA);
               eB = Correlation::et(ks, lengthB, kuhnB);
               correlation[j] += prefactor * x * eA * eB;
            }
         } else {
            for (int j = 0; j < nk; ++j) {
               ks = kSq[j];
               x = std::exp(alpha*ks);
               eA = Correlation::eb(ks, lengthA, kuhnA);
               eB = Correlation::eb(ks, lengthB, kuhnB);
               correlation[j] += prefactor * x * eA * eB;
            }
         }
      }
   }

   /*
   * Increment k-space array of total intramolecular correlation function.
   */
   template <typename WT>
   void Polymer<WT>::computeOmegaTotal(double prefactor, 
                                       Array<double> const & kSq,
                                       Array<double> & correlation) const
   {
      // Preconditions
      UTIL_CHECK(speciesPtr_);
      UTIL_CHECK(nBlock_ > 0);
      UTIL_CHECK(totalLength_ > 0.0);
      int const nk = kSq.capacity();
      UTIL_CHECK(nk > 0);
      UTIL_CHECK(correlation.capacity() == nk);

      // Compute intra-block contributions
      double lengthA, kuhnA, d;
      for (int ia = 0; ia < nBlock_; ++ia) {
         lengthA = length_[ia];
         kuhnA = kuhn_[ia];
         if (PolymerModel::isThread()) {
            for (int j = 0; j < nk; ++j) {
               d = Correlation::dt(kSq[j], lengthA, kuhnA);
               correlation[j] += prefactor * d;
            }
         } else {
            for (int j = 0; j < nk; ++j) {
               d = Correlation::db(kSq[j], lengthA, kuhnA);
               correlation[j] += prefactor * d;
            }
         }
      }

      // Compute inter-block contributions
      if (nBlock_ > 0) {
         double lengthB, kuhnB, alpha, c, ks, x, eA, eB;
         int ia, ib;
         for (ia = 1; ia < nBlock_; ++ia) {
            lengthA = length_[ia];
            kuhnA = kuhn_[ia];
            for (ib = 0; ib < ia; ++ib) {
               lengthB = length_[ib];
               kuhnB = kuhn_[ib];
               alpha = -1.0*rSq_(ia, ib)/6.0;
               c = 2.0 * prefactor;
               if (PolymerModel::isThread()) {
                  for (int j = 0; j < nk; ++j) {
                     ks = kSq[j];
                     x = std::exp(alpha*ks);
                     eA = Correlation::et(ks, lengthA, kuhnA);
                     eB = Correlation::et(ks, lengthB, kuhnB);
                     correlation[j] += c * x * eA * eB;
                  }
               } else {
                  for (int j = 0; j < nk; ++j) {
                     ks = kSq[j];
                     x = std::exp(alpha*ks);
                     eA = Correlation::eb(ks, lengthA, kuhnA);
                     eB = Correlation::eb(ks, lengthB, kuhnB);
                     correlation[j] += c * x * eA * eB;
                  }
               }
            }
         }
      }

   }

} // namespace Correlation
} // namespace Pscf
#endif
