/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HomogeneousComparison.h"

#include <r1d/solvers/Mixture.h>
#include <r1d/solvers/Polymer.h>
#include <r1d/solvers/Block.h>
#include <r1d/solvers/Solvent.h>

#include <pscf/interaction/Interaction.h>
#include <pscf/floryHuggins/FhInteraction.h>
#include <pscf/floryHuggins/FhMixture.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <string>

namespace Pscf {
namespace R1d
{

   using namespace Util;

   /*
   * Default constructor.
   */
   HomogeneousComparison::HomogeneousComparison()
    : SystemAccess(),
      fhMixturePtr_(nullptr),
      fhInteractionPtr_(nullptr),
      hasData_(false)
   {}

   /*
   * Constructor.
   */
   HomogeneousComparison::HomogeneousComparison(System& system)
    : SystemAccess(system)
   {
      fhMixturePtr_ = new FhMixture();
      fhInteractionPtr_ = new FhInteraction();
   }

   /*
   * Destructor.
   */
   HomogeneousComparison::~HomogeneousComparison()
   {
      if (fhMixturePtr_) {
         delete fhMixturePtr_;
      }
      if (fhInteractionPtr_) {
         delete fhInteractionPtr_;
      }
   }

   /*
   * Compute properties of a homogeneous reference system.
   *
   * mode == 0 : Composition equals spatial average composition
   * mode == 1:  Chemical potential equal to that of system,
   *             composition guess given at last grid point.
   * mode == 2:  Chemical potential equal to that of system,
   *             composition guess given at first grid point.
   */
   void HomogeneousComparison::compute(int mode)
   {
    
      UTIL_CHECK(fhMixturePtr_);
      UTIL_CHECK(fhInteractionPtr_);
      fhMixture().initialize(mixture()); 
      fhInteraction() = interaction(); 

      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      int nt = np + ns;
      UTIL_CHECK(fhMixture().nMolecule() == nt);
      if (!p_.isAllocated()) {
         p_.allocate(nt);
      }
      UTIL_CHECK(p_.capacity() == nt);

      if (mode == 0) {

         for (int i = 0; i < np; ++i) {
            p_[i] = mixture().polymer(i).phi();
         }
         fhMixture().setComposition(p_);
         double xi = 0.0;
         fhMixture().computeMu( fhInteraction(), xi );

      } else 
      if (mode == 1 || mode == 2) {

         if (!m_.isAllocated()) {
            m_.allocate(np+ns);
         }
         UTIL_CHECK(m_.capacity() == fhMixture().nMolecule());
         for (int i = 0; i < np; ++i) {
            m_[i] = mixture().polymer(i).mu(); 
         }
         int ix; // Grid index from which we obtain guess of composition
         if (mode == 1) {
            ix = domain().nx() - 1;
         } else 
         if (mode == 2) {
            ix = 0;
         }

         // Compute array p_
         // Define: p_[i] = volume fraction of all blocks of polymer i
         for (int i = 0; i < np; ++i) {
            p_[i] = 0.0;
            int nb = mixture().polymer(i).nBlock();
            for (int j = 0; j < nb; ++j) {
               p_[i] += mixture().polymer(i).block(j).cField()[ix];
            }
         }

         double xi = 0.0;
         fhMixture().computePhi(fhInteraction(), m_, p_, xi);

      } else {
         UTIL_THROW("Unknown mode in computeHomogeneous");
      }

      // Compute Helmholtz free energy and pressure
      fhMixture().computeFreeEnergy(fhInteraction());

      hasData_ = true;
   }

   /*
   * Output properties of homogeneous system and free energy difference.
   *
   * Mode 0:      Outputs fHomo (Helmholtz) and difference df
   * Mode 1 or 2: Outputs Helhmoltz, pressure and difference dp
   */
   void HomogeneousComparison::output(int mode, std::ostream& out) const
   {
      // Precondition
      UTIL_CHECK(hasData_);

      if (mode == 0) {

         // Output free energies
         double fHomo = fhMixture().fHelmholtz();
         double df = system().fHelmholtz() - fHomo;
         out << std::endl;
         out << "f (homo)    = " << Dbl(fHomo, 18, 11) 
                                 << "   [per monomer volume]" << std::endl;
         out << "delta f     = " << Dbl(df, 18, 11)  
                                 << "   [per monomer volume]" << std::endl;

         // Output species-specific properties
         out << std::endl;
         out << "Species:" << std::endl;
         out << "    i"
             << "      mu(homo)      "
             << "      phi(homo)     "
             << std::endl;
         for (int i = 0; i < fhMixture().nMolecule(); ++i) {
            out << Int(i,5)
                << "  " << Dbl(fhMixture().mu(i), 18, 11)
                << "  " << Dbl(fhMixture().phi(i), 18, 11) 
                << std::endl;
         }
         out << std::endl;

      } else
      if (mode == 1 || mode == 2) {

         // Output free energies
         double fHomo = fhMixture().fHelmholtz();
         double pHomo = fhMixture().pressure();
         double pEx   = system().pressure() - pHomo;
         double fEx   = system().fHelmholtz() - fHomo;
         double V     = domain().volume()/mixture().vMonomer();
         double PhiExTot = -1.0*pEx*V;
         double fExTot   = fEx*V; 
         out << std::endl;
         out << "f (homo)   = " << Dbl(fHomo, 18, 11)     
                                << "   [per monomer volume]" << std::endl;
         out << "p (homo)   = " << Dbl(pHomo, 18, 11)     
                                << "   [per monomer volume]" << std::endl;
         out << "delta f    = " << Dbl(fEx, 18, 11)  
                                << "   [per monomer volume]" << std::endl;
         out << "delta p    = " << Dbl(pEx, 18, 11)  
                                << "   [per monomer volume]" << std::endl;
         out << "Phi (ex)   = " << Dbl(PhiExTot, 18, 11) 
                                << "   [total]" << std::endl;
         out << "F (ex)     = " << Dbl(fExTot, 18, 11)  
                                << "   [total]" << std::endl;
         out << "V(tot)/v   = " << Dbl(V, 18, 11) << std::endl;

         // Output species specific properties
         double dV; 
         out << std::endl;
         out << "Species:" << std::endl;
         out << "    i"
             << "      mu            "
             << "      phi(homo)     "
             << "      deltaV        " 
             << std::endl;
         for (int i = 0; i < fhMixture().nMolecule(); ++i) {
            dV = mixture().polymer(i).phi() - fhMixture().phi(i);
            dV *= V;
            out << Int(i,5)
                << "  " << Dbl(fhMixture().mu(i), 18, 11)
                << "  " << Dbl(fhMixture().phi(i), 18, 11) 
                << "  " << Dbl(dV, 18, 11)
                << std::endl;
         }
         out << std::endl;

      }

   }

} // namespace R1d
} // namespace Pscf
