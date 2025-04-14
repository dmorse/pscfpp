#ifndef PRDC_PARAM_ID_CONVERSIONS_H
#define PRDC_PARAM_ID_CONVERSIONS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"

namespace Pscf { 
namespace Prdc { 

   /** 
   * Convert full lattice parameter index to a reduced index.
   * 
   * Here, we define the "full" list of lattice parameters to be
   * {a, b, c, alpha, beta, gamma}, and we define the "reduced"
   * list of lattice parameters to be the list of independent 
   * parameters as defined by a particular lattice system (e.g.,
   * {a, c} for a tetragonal lattice). This function converts an 
   * index of a parameter in the full list, fullId, into the 
   * corresponding index in the reduced list, reducedId. 
   * 
   * If fullId denotes a parameter that is able to vary but is not
   * independent, then the index of the corresponding independent 
   * parameter in the reduced list is returned. For example, in a 
   * tetragonal lattice a == b, so a fullId of 1 would correspond
   * to a reducedId of 0.
   * 
   * If fullId denotes a parameter that is rigid in a given lattice,
   * then reducedId is set to -1. For example, all three  angles are 
   * fixed at 90° in a tetragonal lattice, so a fullId value of 3, 
   * 4, or 5 would result in a reducedId of -1.
   * 
   * \ingroup Prdc_Crystal_Module 
   *
   * \param fullId  the value of fullId to convert
   * \param lattice  the lattice system
   */
   template <int D>
   int convertFullParamIdToReduced(const int fullId, 
                  const typename UnitCell<D>::LatticeSystem lattice);

   /** 
   * Convert reduced lattice parameter index to a full index.
   * 
   * Here, we define the "full" list of lattice parameters to be
   * {a, b, c, alpha, beta, gamma}, and we define the "reduced"
   * list of lattice parameters to be the list of independent 
   * parameters as defined by a particular lattice system (e.g.,
   * {a, c} for a tetragonal lattice). This function converts an 
   * index of a parameter in the reduced list, reducedId, into the 
   * corresponding index in the full list, fullId. 
   * 
   * \ingroup Prdc_Crystal_Module
   *
   * \param reducedId  the value of reducedId to convert
   * \param lattice  the lattice system 
   */
   template <int D>
   int convertReducedParamIdToFull(const int reducedId, 
                  const typename UnitCell<D>::LatticeSystem lattice);


   // Explicit specializations
   template <>
   int convertFullParamIdToReduced<1>(const int fullId, 
                  const typename UnitCell<1>::LatticeSystem lattice);

   template <>
   int convertFullParamIdToReduced<2>(const int fullId, 
                  const typename UnitCell<2>::LatticeSystem lattice);

   template <>
   int convertFullParamIdToReduced<3>(const int fullId, 
                  const typename UnitCell<3>::LatticeSystem lattice);

   template <>
   int convertReducedParamIdToFull<1>(const int reducedId, 
                  const typename UnitCell<1>::LatticeSystem lattice);

   template <>
   int convertReducedParamIdToFull<2>(const int reducedId, 
                  const typename UnitCell<2>::LatticeSystem lattice);

   template <>
   int convertReducedParamIdToFull<3>(const int reducedId, 
                  const typename UnitCell<3>::LatticeSystem lattice);

} // namespace Prdc
} // namespace Pscf
#endif
