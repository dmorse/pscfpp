#ifndef RPC_DOMAIN_H
#define RPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Domain.h>     // base class template

// Forward declarations
namespace Pscf {
   namespace Prdc {
      namespace Cpu {
         template <int D> class WaveList;
         template <int D> class FFT;
      }
   }
   namespace Rpc {
      template <int D> class FieldIo;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Spatial domain for a periodic structure with real fields, on a CPU.
   *
   * Specializations of this template with D =1, 2, and 3 are derived from
   * specializations of the class template Rp::Domain<D, FFT, WLT, FIT>,
   * defined using template type arguments FFT = Prdc::Cpu::FFT\<D\>, 
   * WLT = Prdc::Cpu::WaveList\<D\>, and FIT = Rpc::FieldIo\<D\> that are
   * designed to use standard CPU hardware. The entire public interface 
   * and all of the source code are inherited from this base class. See
   * the documentation of the Rp::Domain base class template for details. 
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Domain 
    : public Rp::Domain< D, FFT<D>, WaveList<D>, FieldIo<D> >
   {};

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations 
namespace Pscf {
   namespace Rp {
      using namespace Prdc;
      extern template 
      class Domain<1, Cpu::FFT<1>, Cpu::WaveList<1>, Rpc::FieldIo<1> >;
      extern template 
      class Domain<2, Cpu::FFT<2>, Cpu::WaveList<2>, Rpc::FieldIo<2> >;
      extern template 
      class Domain<3, Cpu::FFT<3>, Cpu::WaveList<3>, Rpc::FieldIo<3> >;
   } 
   namespace Rpc {
      extern template class Domain<1>;
      extern template class Domain<2>;
      extern template class Domain<3>;
   } 
} 
#endif
