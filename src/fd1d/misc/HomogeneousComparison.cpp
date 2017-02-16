/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HomogeneousComparison.h"

#include <pscf/inter/Interaction.h>
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <string>
//#include <ctime>
//#include <iomanip>
//#include <sstream>
//#include <unistd.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /*
   * Default constructor.
   */
   HomogeneousComparison::HomogeneousComparison()
    : SystemAccess()
   {}

   /*
   * Constructor.
   */
   HomogeneousComparison::HomogeneousComparison(System& system)
    : SystemAccess(system)
   {}

   /*
   * Destructor.
   */
   HomogeneousComparison::~HomogeneousComparison()
   {}

   #if 0
   void HomogeneousComparison::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      //int ns = mixture().nSolvent(); 
      int ns = 0;
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      if (c_.isAllocated()) {
         UTIL_CHECK(c_.capacity() == nm);
      } else {
         c_.allocate(nm);
      }
      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      int i;   // molecule index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      for (i = 0; i < np; ++i) {

         // Initial array of clump sizes 
         for (j = 0; j < nm; ++j) {
            c_[j] = 0.0;
         }

         // Compute clump sizes for all monomer types.
         nb = mixture().polymer(i).nBlock(); 
         for (k = 0; k < nb; ++k) {
            Block& block = mixture().polymer(i).block(k);
            j = block.monomerId();
            c_[j] += block.length();
         }
 
         // Count the number of clumps of nonzero size
         nc = 0;
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 1.0E-8) {
               ++nc;
            }
         }
         homogeneous_.molecule(i).setNClump(nc);
 
         // Set clump properties for this Homogeneous::Molecule
         k = 0; // Clump index
         for (j = 0; j < nm; ++j) {
            if (c_[j] > 1.0E-8) {
               homogeneous_.molecule(i).clump(k).setMonomerId(j);
               homogeneous_.molecule(i).clump(k).setSize(c_[j]);
               ++k;
            }
         }
         homogeneous_.molecule(i).computeSize();

         #if 0
         {
            std::cout << "Molecule # " << i << std::endl;
            nc = homogeneous_.molecule(i).nClump();
            std::cout << "nClump = " << nc << std::endl;
            double size;
            for (k = 0; k < nc; ++k) {
               j = homogeneous_.molecule(i).clump(k).monomerId();
               size = homogeneous_.molecule(i).clump(k).size();
               std::cout << k << "  " << j << "  " << size << "\n";
            }
         }
         #endif
        
      }

   }
   #endif

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
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      if (!p_.isAllocated()) {
         p_.allocate(np+ns);
      }
      UTIL_CHECK(p_.capacity() == homogeneous().nMolecule());

      if (mode == 0) {

         for (int i = 0; i < np; ++i) {
            p_[i] = mixture().polymer(i).phi();
         }
         homogeneous().setComposition(p_);
         double xi = 0.0;
         homogeneous().computeMu(interaction(), xi);

      } else 
      if (mode == 1 || mode == 2) {

         if (!m_.isAllocated()) {
            m_.allocate(np+ns);
         }
         UTIL_CHECK(m_.capacity() == homogeneous().nMolecule());
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

         for (int i = 0; i < np; ++i) {
            p_[i] = 0.0;
            int nb = mixture().polymer(i).nBlock();
            for (int j = 0; j < nb; ++j) {
               p_[i] += mixture().polymer(i).block(j).cField()[ix];
            }
         }

         #if 0
         std::cout << std::endl;
         std::cout << "Composition at boundary " << std::endl;
         for (int i = 0; i < np; ++i) {
            std::cout << "phi[" << i << "] = " << p_[i] << std::endl;
         }
         #endif
    
         double xi = 0.0;
         homogeneous().computePhi(interaction(), m_, p_, xi);

      } else {
         UTIL_THROW("Unknown mode in computeHomogeneous");
      }

      // Compute Helmholtz free energy and pressure
      homogeneous().computeFreeEnergy(interaction());
   }

   /*
   * Output properties of homogeneous system and free energy difference.
   *
   * Mode 0:      Outputs fHomo (Helmholtz) and difference df
   * Mode 1 or 2: Outputs Helhmoltz, pressure and difference dp
   */
   void HomogeneousComparison::output(int mode, std::ostream& out)
   {
      // Output free energies
      out << std::endl;
      if (mode == 0) {

         // Output free energies
         double fHomo = homogeneous().fHelmholtz();
         double df = system().fHelmholtz() - fHomo;
         out << "f (homo)    = " << Dbl(fHomo, 18, 11) << std::endl;
         out << "delta f     = " << Dbl(df, 18, 11)    << std::endl;

         // Output polymer properties
         out << std::endl;
         out << "Polymers: i, mu(homo)[i], phi[i] " << std::endl;
         for (int i = 0; i < homogeneous().nMolecule(); ++i) {
            out << i  
                << "  " << Dbl(homogeneous().mu(i), 18, 11)
                << "  " << Dbl(homogeneous().phi(i), 18, 11) 
                << std::endl;
         }
         out << std::endl;

      } else
      if (mode == 1 || mode == 2) {

         // Output free energies
         double fHomo = homogeneous().fHelmholtz();
         double pHomo = homogeneous().pressure();
         double pEx   = system().pressure() - pHomo;
         double fEx   = system().fHelmholtz() - fHomo;
         double V     = domain().volume()/mixture().vMonomer();
         double PExV  = -1.0*pEx*V;
         double FExV  = fEx*V; 
         out << "f (homo)   = " << Dbl(fHomo, 18, 11) << std::endl;
         out << "p (homo)   = " << Dbl(pHomo, 18, 11) << std::endl;
         out << "f (ex)     = " << Dbl(fEx, 18, 11)   << std::endl;
         out << "p (homo)   = " << Dbl(pEx, 18, 11)   << std::endl;
         out << "-p(ex)*V   = " << Dbl(PExV, 18, 11)  << std::endl;
         out << "f(ex)*V    = " << Dbl(FExV, 18, 11)  << std::endl;

         // Output polymer properties
         double dV; 
         out << std::endl;
         out << "Polymers:" << std::endl;
         out << "    i"
             << "      mu[i]         "
             << "      phi(homo)[i]  "
             << "      dV[i]         " 
             << std::endl;
         for (int i = 0; i < homogeneous().nMolecule(); ++i) {
            dV = mixture().polymer(i).phi() - homogeneous().phi(i);
            dV *= V;
            out << Int(i,5)
                << "  " << Dbl(homogeneous().mu(i), 18, 11)
                << "  " << Dbl(homogeneous().phi(i), 18, 11) 
                << "  " << Dbl(dV, 18, 11)
                << std::endl;
         }
         out << std::endl;

      }

   }

} // namespace Fd1d
} // namespace Pscf
