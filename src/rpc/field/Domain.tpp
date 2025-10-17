#ifndef RPC_DOMAIN_TPP
#define RPC_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <prdc/field/DomainTmpl.tpp>

#include <rpc/field/FieldIo.tpp>
#include <prdc/cpu/WaveList.h>
#include <prdc/cpu/FFT.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : DomainTmpl<D, FFT<D>, WaveList<D>, FieldIo<D> >()
   {  ParamComposite::setClassName("Domain"); }

} // namespace Rpc
} // namespace Pscf
#endif
