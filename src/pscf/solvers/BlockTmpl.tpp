#ifndef PSCF_BLOCK_TMPL_TPP
#define PSCF_BLOCK_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockTmpl.h"

namespace Pscf {

   /*
   * Constructor.
   */
   template <class QT, class FT>
   BlockTmpl<QT,FT>::BlockTmpl()
    : propagators_(),
      cField_(),
      kuhn_(0.0)
   {
      propagators_.allocate(2);
      propagator(0).setDirectionId(0);
      propagator(1).setDirectionId(1);
      propagator(0).setPartner(propagator(1));
      propagator(1).setPartner(propagator(0));
   }

   /*
   * Destructor.
   */
   template <class QT, class FT>
   BlockTmpl<QT,FT>::~BlockTmpl()
   {}

   /*
   * Set the monomer statistical segment length.
   */
   template <class QT, class FT>
   void BlockTmpl<QT,FT>::setKuhn(double kuhn)
   {  kuhn_ = kuhn; }

}
#endif
