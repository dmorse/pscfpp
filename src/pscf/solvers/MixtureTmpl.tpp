#ifndef PSCF_MIXTURE_TMPL_TPP
#define PSCF_MIXTURE_TMPL_TPP

/*
* PSCF - Mixture Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureTmpl.h"

namespace Pscf
{

   /*
   * Constructor.
   */
   template <class PT, class ST>
   MixtureTmpl<PT,ST>::MixtureTmpl()
    : MixtureBase(),
      ParamComposite(),
      polymers_(),
      solvents_()
   {}

   /*
   * Destructor.
   */
   template <class PT, class ST>
   MixtureTmpl<PT,ST>::~MixtureTmpl()
   {}

   /*
   * Get a PolymerSpecies descriptor by non-const reference.
   */
   template <class PT, class ST>
   PolymerSpecies const & MixtureTmpl<PT,ST>::polymerSpecies(int id) const
   {
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   /*
   * Get a SolventSpecies descriptor by const reference.
   */
   template <class PT, class ST>
   SolventSpecies const & MixtureTmpl<PT,ST>::solventSpecies(int id) const
   {
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id];
   }

   /*
   * Read all parameters and initialize.
   */
   template <class PT, class ST>
   void MixtureTmpl<PT,ST>::readParameters(std::istream& in)
   {
      // Read nMonomer and monomers array
      read<int>(in, "nMonomer", nMonomer_);
      monomers_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         monomers_[i].setId(i);
      }
      readDArray< Monomer >(in, "monomers", monomers_, nMonomer_);

      /*
      * The input format for a single monomer is defined in the istream
      * extraction operation (operator >>) for a Pscf::Monomer, in file
      * pscf/chem/Monomer.cpp. The text representation contains only the
      * value for the monomer statistical segment Monomer::kuhn.
      */

      // Read nPolymer (required parameter, must be > 0)
      read<int>(in, "nPolymer", nPolymer_);
      UTIL_CHECK(nPolymer_ > 0);

      // Optionally read nSolvent, with nSolvent=0 by default
      nSolvent_ = 0;
      readOptional<int>(in, "nSolvent", nSolvent_);

      // Read polymers and compute nBlock
      nBlock_ = 0;
      if (nPolymer_ > 0) {

         polymers_.allocate(nPolymer_);
         for (int i = 0; i < nPolymer_; ++i) {
            readParamComposite(in, polymer(i));
            nBlock_ = nBlock_ + polymer(i).nBlock();
         }

         // Set statistical segment lengths for all blocks
         double kuhn;
         int monomerId;
         for (int i = 0; i < nPolymer_; ++i) {
            for (int j = 0; j < polymer(i).nBlock(); ++j) {
               monomerId = polymer(i).block(j).monomerId();
               kuhn = monomer(monomerId).kuhn();
               polymer(i).block(j).setKuhn(kuhn);
            }
         }

      }

      // Read solvents
      if (nSolvent_ > 0) {

         solvents_.allocate(nSolvent_);
         for (int i = 0; i < nSolvent_; ++i) {
            readParamComposite(in, solvent(i));
         }

      }

      // Optionally read monomer reference value
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);

   }

}
#endif
