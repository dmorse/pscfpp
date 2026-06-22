#ifndef RP_MIXTURE_MODIFIER_TPP
#define RP_MIXTURE_MODIFIER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureModifier.h"
#include <pscf/chem/Monomer.h>

namespace Pscf {
namespace Rp {

   /*
   * Constructor
   */
   template <int D, class T>
   MixtureModifier<D,T>::MixtureModifier()
    : mixturePtr_(nullptr)
   {}

   /*
   * Destructor
   */
   template <int D, class T>
   MixtureModifier<D,T>::~MixtureModifier()
   {}

   /*
   * Create an association with a mixture.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::associate(Mixture<D,T>& mixture)
   {
      UTIL_CHECK(!mixturePtr_);  
      mixturePtr_ = &mixture; 
   }

   // Parameter Modifiers

   /*
   * Set statistical segment length for one monomer type.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setKuhn(int monomerId, double kuhn)
   {  mixture().setKuhn(monomerId, kuhn); }

   /*
   * Set volume fraction for a polymer.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setPhiPolymer(int polymerId, 
                                               double phi)
   {  mixture().polymer(polymerId).setPhi(phi); }

   /*
   * Set chemical potential for a polymer.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setMuPolymer(int polymerId, 
                                               double mu)
   {  mixture().polymer(polymerId).setMu(mu); }

   /*
   * Set the length of a polymer block.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setBlockLength(int polymerId, 
                                                int blockId,
                                                double length)
   {  mixture().polymer(polymerId).block(blockId).setLength(length); }

   /*
   * Set the volume fraction for a solvent.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setPhiSolvent(int solventId, 
                                               double phi)
   {  mixture().solvent(solventId).setPhi(phi); }

   /*
   * Set the chemical potential for a solvent species.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setMuSolvent(int solventId, 
                                              double mu)
   {  mixture().solvent(solventId).setMu(mu); }

   /*
   * Set the size of solvent species.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setSolventSize(int solventId, 
                                                double size)
   {  mixture().solvent(solventId).setSize(size); }

   /*
   * Set the monomer reference volume for the mixture.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::setVMonomer(double vMonomer)
   {  mixture().setVmonomer(vMonomer); }

   // Other public non-const functions

   /*
   * Clear all data that depends on the unit cell parameters.
   */
   template <int D, class T>
   void MixtureModifier<D,T>::clearUnitCellData()
   {  mixture().clearUnitCellData(); }

   // Private memmber function

   /*
   * Get the associated mixture by reference
   */
   template <int D, class T>
   Mixture<D,T>& MixtureModifier<D,T>::mixture()
   {
      UTIL_CHECK(mixturePtr_);
      return *mixturePtr_;  
   }

} // namespace Rp
} // namespace Pscf
#endif
