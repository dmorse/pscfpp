/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <fd1d/domain/Domain.h>

#include <cmath>

namespace Pscf { 
namespace Fd1d
{ 

   Mixture::Mixture()
    : ds_(-1.0),
      domainPtr_(0)
   {  setClassName("Mixture"); }

   Mixture::~Mixture()
   {}

   void Mixture::readParameters(std::istream& in)
   {
      MixtureTmpl<Polymer, Solvent>::readParameters(in);

      // Read optimal contour step size
      read(in, "ds", ds_);

      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
   }

   void Mixture::setDomain(Domain const& domain)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);

      domainPtr_ = &domain;

      // Process polymers - set discretization for all blocks
      if (nPolymer() > 0) {
         int i, j;
         for (i = 0; i < nPolymer(); ++i) {
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               polymer(i).block(j).setDiscretization(domain, ds_);
            }
         }
      }

      // Process solvents - set discretization for all solvents
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).setDiscretization(domain);
         }
      }

   }

   /*
   * Reset statistical segment length for one monomer type.
   */
   void Mixture::setKuhn(int monomerId, double kuhn)
   {
      // Set new Kuhn length for relevant Monomer object
      monomer(monomerId).setKuhn(kuhn);

      // Update kuhn length for all blocks of this monomer type
      for (int i = 0; i < nPolymer(); ++i) {
         for (int j =  0; j < polymer(i).nBlock(); ++j) {
            if (monomerId == polymer(i).block(j).monomerId()) {
               polymer(i).block(j).setKuhn(kuhn);
            }
         }
      }
   }

   /*
   * Compute concentrations (but not total free energy).
   */
   void Mixture::compute(DArray<Mixture::WField> const & wFields, 
                         DArray<Mixture::CField>& cFields)
   {
      UTIL_CHECK(domainPtr_);
      UTIL_CHECK(domain().nx() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

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

      // Process polymer species
      if (nPolymer() > 0) {

         for (i = 0; i < nPolymer(); ++i) {

            // Solve MDE for all blocks in polymer
            polymer(i).compute(wFields);

            // Accumulate monomer concentrations
            for (j = 0; j < polymer(i).nBlock(); ++j) {
               int monomerId = polymer(i).block(j).monomerId();
               UTIL_CHECK(monomerId >= 0);
               UTIL_CHECK(monomerId < nm);
               CField& monomerField = cFields[monomerId];
               CField const & blockField = polymer(i).block(j).cField();
               for (k = 0; k < nx; ++k) {
                  monomerField[k] += blockField[k];
               }
            }

         }

      }

      // Process solvent species
      if (nSolvent() > 0) {
         for (i = 0; i < nSolvent(); ++i) {

            int monomerId = solvent(i).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);

            // Compute solvent concentration and q
            solvent(i).compute(wFields[monomerId]);

            // Add to monomer concentrations
            CField& monomerField = cFields[monomerId];
            CField const & solventField = solvent(i).cField();
            for (k = 0; k < nx; ++k) {
               monomerField[k] += solventField[k];
            }

         }
      }

   }

} // namespace Fd1d
} // namespace Pscf
