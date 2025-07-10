#ifndef PRDC_SCFT_REAL_TPP
#define PRDC_SCFT_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ScftReal.h"

#include <pscf/inter/Interaction.h>
#include <pscf/mesh/Mesh.h>

//#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

//#include <pscf/math/IntVec.h>
//#include <prdc/crystal/UnitCell.h>
//#include <util/containers/DArray.h>
//#include <util/containers/FSArray.h>
//#include <util/misc/ioUtil.h>
//#include <string>

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D, class ST>
   ScftReal<D,ST>::ScftReal(SystemT const & system)
    : SystemConstRefReal<ST>(system),
      fHelmholtz_(0.0),
      fIdeal_(0.0),
      fInter_(0.0),
      fExt_(0.0),
      pressure_(0.0),
      hasFreeEnergy_(false)
   {}

   /*
   * Destructor.
   */
   template <int D, class ST>
   ScftReal<D,ST>::~ScftReal()
   {}

   /*
   * Compute Helmholtz free energy and pressure.
   */
   template <int D, class ST>
   void ScftReal<D,ST>::compute()
   {
      if (hasFreeEnergy_) return;

      UTIL_CHECK(w().hasData());
      UTIL_CHECK(c().hasData());

      int nm = mixture().nMonomer();   // number of monomer types
      int np = mixture().nPolymer();   // number of polymer species
      int ns = mixture().nSolvent();   // number of solvent species

      // Initialize all free energy contributions to zero
      fHelmholtz_ = 0.0;
      fIdeal_ = 0.0;
      fInter_ = 0.0;
      fExt_ = 0.0;

      double phi, mu;

      // Compute polymer ideal gas contributions to fIdeal_
      if (np > 0) {
         PolymerSpecies const * polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture().polymerSpecies(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            if (PolymerModel::isThread()) {
               length = polymerPtr->length();
            } else {
               length = (double) polymerPtr->nBead();
            }
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 ) / length;
               // Recall: mu = ln(phi/q)
            }
         }
      }

      // Compute solvent ideal gas contributions to fIdeal_
      if (ns > 0) {
         SolventSpecies const * solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture().solventSpecies(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 )/size;
            }
         }
      }

      // Volume integrals with a mask: If the system has a mask, then the
      // volume that should be used in calculating free energy/pressure
      // is the volume available to the material, not the total unit cell
      // volume. We thus divide all terms that involve integrating over
      // the unit cell volume by quantity mask().phiTot(), which is the
      // volume fraction of the unit cell that is occupied by material.
      // This properly scales them to the correct value. fExt_, fInter_,
      // and the Legendre transform component of fIdeal_ all require
      // this scaling. If no mask is present, mask().phiTot() = 1 and no
      // scaling occurs.
      const double phiTot = mask().phiTot();

      // Compute Legendre transform subtraction from fIdeal_
      double temp = 0.0;
      if (w().isSymmetric()) {
         // Use expansion in symmetry-adapted orthonormal basis
         UTIL_CHECK(w().isAllocatedBasis());
         UTIL_CHECK(c().isAllocatedBasis());
         const int nBasis = domain().basis().nBasis();
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               temp -= w().basis(i)[k] * c().basis(i)[k];
            }
         }
      } else {
         // Use summation over grid points
         const int meshSize = domain().mesh().size();
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < meshSize; ++k) {
               temp -= w().rgrid(i)[k] * c().rgrid(i)[k];
            }
         }
         temp /= double(meshSize);
      }
      temp /= phiTot;
      fIdeal_ += temp;
      fHelmholtz_ += fIdeal_;

      // Compute contribution from external fields, if they exist
      if (system().hasExternalFields()) {
         if (w().isSymmetric() && h().isSymmetric()) {
            // Use expansion in symmetry-adapted orthonormal basis
            UTIL_CHECK(h().isAllocatedBasis());
            UTIL_CHECK(c().isAllocatedBasis());
            const int nBasis = domain().basis().nBasis();
            for (int i = 0; i < nm; ++i) {
               for (int k = 0; k < nBasis; ++k) {
                  fExt_ += h().basis(i)[k] * c().basis(i)[k];
               }
            }
         } else {
            // Use summation over grid points
            const int meshSize = domain().mesh().size();
            for (int i = 0; i < nm; ++i) {
               for (int k = 0; k < meshSize; ++k) {
                  fExt_ += h().rgrid(i)[k] * c().rgrid(i)[k];
               }
            }
            fExt_ /= double(meshSize);
         }
         fExt_ /= phiTot;
         fHelmholtz_ += fExt_;
      }

      // Compute excess interaction free energy [ phi^{T}*chi*phi/2 ]
      if (w().isSymmetric()) {
         UTIL_CHECK(c().isAllocatedBasis());
         const int nBasis = domain().basis().nBasis();
         for (int i = 0; i < nm; ++i) {
            for (int j = i; j < nm; ++j) {
               const double chi = interaction().chi(i,j);
               if (std::abs(chi) > 1.0E-9) {
                  double temp = 0.0;
                  for (int k = 0; k < nBasis; ++k) {
                     temp += c().basis(i)[k] * c().basis(j)[k];
                  }
                  if (i == j) {
                     fInter_ += 0.5*chi*temp;
                  } else {
                     fInter_ += chi*temp;
                  }
               }
            }
         }
      } else {
         const int meshSize = domain().mesh().size();
         for (int i = 0; i < nm; ++i) {
            for (int j = i; j < nm; ++j) {
               const double chi = interaction().chi(i,j);
               if (std::abs(chi) > 1.0E-9) {
                  double temp = 0.0;
                  for (int k = 0; k < meshSize; ++k) {
                     temp += c().rgrid(i)[k] * c().rgrid(j)[k];
                  }
                  if (i == j) {
                     fInter_ += 0.5*chi*temp;
                  } else {
                     fInter_ += chi*temp;
                  }
               }
            }
         }
         fInter_ /= double(meshSize);
      }
      fInter_ /= phiTot;
      fHelmholtz_ += fInter_;

      // Initialize pressure (-1 x grand-canonical free energy / monomer)
      pressure_ = -fHelmholtz_;

      // Polymer chemical potential corrections to pressure
      if (np > 0) {
         PolymerSpecies const * polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture().polymerSpecies(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            if (PolymerModel::isThread()) {
               length = polymerPtr->length();
            } else {
               length = (double) polymerPtr->nBead();
            }
            if (phi > 1.0E-08) {
               pressure_ += mu * phi / length;
            }
         }
      }

      // Solvent corrections to pressure
      if (ns > 0) {
         SolventSpecies const * solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture().solventSpecies(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1.0E-08) {
               pressure_ += mu * phi / size;
            }
         }
      }

      hasFreeEnergy_ = true;
   }

   /*
   * Clear stored values of thermodynamic properties.
   */
   template <int D, class ST>
   void ScftReal<D, ST>::clear()
   {  hasFreeEnergy_ = false; }

   /*
   * Write thermodynamic properties to file.
   */
   template <int D, class ST>
   void ScftReal<D, ST>::write(std::ostream& out)
   {
      if (!hasFreeEnergy_) {
         compute();
      }

      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;
      out << "fIdeal        " << Dbl(fIdeal_, 18, 11) << std::endl;
      out << "fInter        " << Dbl(fInter_, 18, 11) << std::endl;
      if (system().hasExternalFields()) {
         out << "fExt          " << Dbl(fExt_, 18, 11) << std::endl;
      }
      out << std::endl;

      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();

      if (np > 0) {
         out << "polymers:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          "
             << std::endl;
         for (int i = 0; i < np; ++i) {
            out << Int(i, 5)
                << "  " << Dbl(mixture().polymer(i).phi(),18, 11)
                << "  " << Dbl(mixture().polymer(i).mu(), 18, 11)
                << std::endl;
         }
         out << std::endl;
      }

      if (ns > 0) {
         out << "solvents:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          "
             << std::endl;
         for (int i = 0; i < ns; ++i) {
            out << Int(i, 5)
                << "  " << Dbl(mixture().solvent(i).phi(),18, 11)
                << "  " << Dbl(mixture().solvent(i).mu(), 18, 11)
                << std::endl;
         }
         out << std::endl;
      }

      out << "cellParams:" << std::endl;
      for (int i = 0; i < domain().unitCell().nParameter(); ++i) {
         out << Int(i, 5)
             << "  "
             << Dbl(domain().unitCell().parameter(i), 18, 11)
             << std::endl;
      }
      out << std::endl;
   }

} // namespace Prdc
} // namespace Pscf
#endif
