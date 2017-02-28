#ifndef PSSP_MIXTURE_TPP
#define PSSP_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <fd1d/domain/Domain.h>

#include <cmath>

namespace Pscf { 
namespace Fd1d
{ 

   template <int D>
   Mixture<D>::Mixture()
    : vMonomer_(1.0),
      ds_(-1.0),
      domainPtr_(0)
   {  setClassName("Mixture"); }

   template <int D>
   Mixture<D>::~Mixture()
   {}

   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureTmpl<Polymer, Solvent>::readParameters(in);
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);
      read(in, "ds", ds_);

      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
   }

   template <int D>
   void Mixture<D>::setDomain(Domain const& domain)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);

      domainPtr_ = &domain;

      // Set discretization for all blocks
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).setDiscretization(domain, ds_);
         }
      }

   }

   /*
   * Compute concentrations (but not total free energy).
   */
   template <int D>
   void Mixture<D>::compute(DArray<Mixture<D>::WField> const & wFields, 
                            DArray<Mixture<D>::CField>& cFields)
   {
      //UTIL_CHECK(domainPtr_);
      //UTIL_CHECK(domain().nx() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

      #if 0
      int nx = domain().nx();
      int nm = nMonomer();
      int i, j, k;

      // Clear all monomer concentration fields
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == nx);
         UTIL_CHECK(wFields[i].capacity() == nx);
         for (j = 0; j < nx; ++j) {
            cFields[i][j] = 0.0;
         }
      }

      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields);
      }

      // Accumulate monomer concentration fields
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            int monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            CField& monomerField = cFields[monomerId];
            CField& blockField = polymer(i).block(j).cField();
            for (k = 0; k < nx; ++k) {
               monomerField[k] += blockField[k];
            }
         }
      }
      #endif

      // To do: Add compute functions and accumulation for solvents.

   }

} // namespace Fd1d
} // namespace Pscf
#endif
