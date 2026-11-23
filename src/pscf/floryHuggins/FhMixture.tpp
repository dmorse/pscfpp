#ifndef PSCF_FLORY_HUGGINS_MIXTURE_TPP
#define PSCF_FLORY_HUGGINS_MIXTURE_TPP

/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
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

   /*
   * Initialize state to that of a Composition<W>.
   *
   * This function allocates iff necessary, and may be called more
   * than once.
   */
   template <typename WT>
   void FhMixture::initialize(Composition<WT> const & mixture)
   {

      // Set number of molecular species and monomers
      int nm = mixture.nMonomer(); 
      int np = mixture.nPolymer(); 
      int ns = mixture.nSolvent(); 
      int nt = np + ns;

      // Check nMonomer, set if necessary
      if (nMonomer_ == 0) {
         setNMonomer(nm);
      } 
      UTIL_CHECK(nMonomer() == nm);
      UTIL_CHECK(nMonomer_ == nm);
      UTIL_CHECK(c_.capacity() == nm);
      UTIL_CHECK(w_.capacity() == nm);

      // Check nMolecule, set if necessary
      if (nMolecule_ == 0) {
         setNMolecule(nt);
      }
      UTIL_CHECK(nMolecule_ == nt);
      UTIL_CHECK(molecules_.capacity() == nt);
      UTIL_CHECK(phi_.capacity() == nt);
      UTIL_CHECK(mu_.capacity() == nt);

      /*
      * Work space for identifying clumps for each polymer species
      * clumpSizes[j] is the sum of block lengths of monomer type j
      * clumpFlags[j] is true iff a molecule contain a block of type j
      */
      DArray<double> clumpSizes;
      DArray<bool> clumpFlags;
      clumpSizes.allocate(nm);
      clumpFlags.allocate(nm);

      int i;   // species index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks in a polymer species
      int nc;  // number of clumps in a polymer species
 
      // Loop over polymer molecule species
      if (np > 0) {
         for (i = 0; i < np; ++i) {
             PolymerSpecies<WT> const& polymer = mixture.polymerSpecies(i);
   
            // Initial array of clump sizes 
            for (j = 0; j < nm; ++j) {
               clumpSizes[j] = 0.0;
               clumpFlags[j] = false;
            }
 
            // Compute clump sizes and flags for all monomer types.
            nb = polymer.nBlock(); 
            for (k = 0; k < nb; ++k) {
               Edge const& edge = polymer.edge(k);
               UTIL_CHECK(edge.length() > 0);
               j = edge.monomerId();
               clumpSizes[j] += edge.length();
               clumpFlags[j] = true;
            }
    
            // Count the number nc of clumps
            nc = 0;
            for (j = 0; j < nm; ++j) {
               UTIL_CHECK(clumpSizes[j] >= 0.0);
               if (clumpFlags[j]) {
                  UTIL_CHECK(clumpSizes[j] > 0.0);
                  ++nc;
               } else {
                  UTIL_CHECK(clumpSizes[j] == 0.0);
               }
            }
            molecule(i).setNClump(nc);
    
            // Set clump properties for this FhMolecule species
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (clumpFlags[j]) {
                  molecule(i).clump(k).setMonomerId(j);
                  molecule(i).clump(k).setSize(clumpSizes[j]);
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
