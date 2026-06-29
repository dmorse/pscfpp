/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"               // Class header
#include <cpc/field/FieldIo.tpp>
#include <prdc/field/cpu/WaveList.h>
#include <prdc/field/cpu/FFT.h>

#include <cp/field/Domain.tpp>    // Base class template implementation

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : Cp::Domain<D, FFT<D>, WaveList<D>, FieldIo<D> >()
   {  ParamComposite::setClassName("Domain"); }

} // namespace Cpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf {
   namespace Cp {
      using namespace Prdc::Cpu;
      template class Domain<1, FFT<1>, WaveList<1>, Cpc::FieldIo<1> >;
      template class Domain<2, FFT<2>, WaveList<2>, Cpc::FieldIo<2> >;
      template class Domain<3, FFT<3>, WaveList<3>, Cpc::FieldIo<3> >;
   } 
   namespace Cpc {
      template class Domain<1>;
      template class Domain<2>;
      template class Domain<3>;
   } 
} 
