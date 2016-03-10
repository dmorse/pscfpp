/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include "Domain.h"

#include <cmath>

namespace Pscf { 
namespace Fd1d
{ 

   Mixture::Mixture()
   {  setClassName("Mixture"); }

   Mixture::~Mixture()
   {}

   void Mixture::readParameters(std::istream& in)
   {
      MixtureTmpl<Polymer, Solvent>::readParameters(in);
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

      // Set discretization for all blocks
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).setDiscretization(domain, ds_);
         }
      }

   }

   void Mixture::compute(DArray<Mixture::WField> const & wFields, 
                         DArray<Mixture::CField>& cFields, bool thermoFlag)
   {
      UTIL_CHECK(domainPtr_);
      UTIL_CHECK(domain().nx() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

      int i, j, k;

      // Clear all monomer concentration fields
      int nx = domain().nx();
      for (i = 0; i < nMonomer(); ++i) {
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
            UTIL_CHECK(monomerId < nMonomer());
            CField& monomerField = cFields[monomerId];
            CField& blockField = polymer(i).block(j).cField();
            for (k = 0; k < nx; ++k) {
               monomerField[k] += blockField[k];
            }
         }
      }

      if (thermoFlag) {
         fHelmholtz_ = 0.0;
         double phi, mu, length;
         for (i = 0; i < nPolymer(); ++i) {
            phi = polymer(i).phi();
            mu = polymer(i).mu();
            length = polymer(i).length();
            fHelmholtz_ += phi*( log(mu) - 1.0 )/length;
         }
         for (i = 0; i < nMonomer(); ++i) {
            fHelmholtz_ -= domain().innerProduct(wFields[i],cFields[i]);
         }
         if (!work_.isAllocated()) {
            work_.allocate(nx);
         }
         #if 0
         for (i = 0; i < nx; ++i) {
            for (j = 0; i < nMonomer(); ++i) {
            }
            work_[i] = interaction().fHelmholtz(cArray);
         }
         #endif
      }
    
   }

}
}
