#ifndef PSCF_FLORY_HUGGINS_MIXTURE_TPP
#define PSCF_FLORY_HUGGINS_MIXTURE_TPP

/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FhMixture.h"
#include "FhMolecule.h"
#include "FhClump.h"
#include <pscf/chem/Composition.h>
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/SolventSpecies.h>
#include <pscf/chem/Edge.h>

namespace Pscf {

   using namespace Util;

   template <typename WT>
   void FhMixture::initialize(Composition<WT> const & mixture)
   {

      // Set number of molecular species and monomers
      int nm = mixture.nMonomer(); 
      int np = mixture.nPolymer(); 
      int ns = mixture.nSolvent(); 
      int nt = np + ns;

      // Check nMonomer, set if necessary
      if (nMonomer() == 0) {
         setNMonomer(nm);
      } 
      UTIL_CHECK(nMonomer() == nm);
      UTIL_CHECK(c_.capacity() == nm);
      UTIL_CHECK(w_.capacity() == nm);

      // Check nMolecule, set if necessary
      if (nMolecule() == 0) {
         setNMolecule(nt);
      }
      UTIL_CHECK(nMolecule() == nt);
      UTIL_CHECK(molecules_.capacity() == nt);
      UTIL_CHECK(phi_.capacity() == nt);
      UTIL_CHECK(mu_.capacity() == nt);

      // Work space of clump sizes
      DArray<double> c_;
      c_.allocate(nm);

      int i;   // species index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      if (np > 0) {
         for (i = 0; i < np; ++i) {
             PolymerSpecies<WT> const& polymer = mixture.polymerSpecies(i);
   
            // Initial array of clump sizes 
            for (j = 0; j < nm; ++j) {
               c_[j] = 0.0;
            }
   
            // Compute clump sizes for all monomer types.
            nb = polymer.nBlock(); 
            for (k = 0; k < nb; ++k) {
               Edge const& edge = polymer.edge(k);
               j = edge.monomerId();
               c_[j] += edge.length();
            }
    
            // Count the number of clumps of nonzero size
            nc = 0;
            for (j = 0; j < nm; ++j) {
               if (c_[j] > 1.0E-8) {
                  ++nc;
               }
            }
            molecule(i).setNClump(nc);
    
            // Set clump properties for this FhMolecule
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (c_[j] > 1.0E-8) {
                  molecule(i).clump(k).setMonomerId(j);
                  molecule(i).clump(k).setSize(c_[j]);
                  ++k;
               }
            }
            UTIL_CHECK(k == nc);
            molecule(i).computeSize();
   
         }
      }

      // Add solvent contributions
      if (ns > 0) {
         double size;
         int monomerId;
         for (int is = 0; is < ns; ++is) {
            i = is + np;
            SolventSpecies<WT> const& solvent = mixture.solventSpecies(is);
            monomerId = solvent.monomerId();
            size = solvent.size();
            molecule(i).setNClump(1);
            molecule(i).clump(0).setMonomerId(monomerId);
            molecule(i).clump(0).setSize(size);
            molecule(i).computeSize();
         }
      }

   }

} // namespace Pscf
#endif
