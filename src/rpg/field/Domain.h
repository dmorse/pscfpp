#ifndef RPG_DOMAIN_H
#define RPG_DOMAIN_H

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
      namespace Cuda {
         template <int D> class FFT;
         template <int D> class WaveList;
      }
   }
   namespace Rpg {
      template <int D> class FieldIo;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Spatial domain for a periodic structure with real fields, on a GPU.
   *
   * See the interface of the Rp::Domain base class template for
   * complete API documentation. The Rpg::Domain class template is 
   * simply a named partial specialization of the base class template, 
   * defined using template type parameters FFT = Prdc::Cuda::FFT<D>, 
   * WLT = Prdc::Cuda::WaveList<D>, and FIT = Rpg::FieldIo<D> . The 
   * public interface is identical to that of the base class.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Domain 
     : public Rp::Domain< D, FFT<D>, WaveList<D>, FieldIo<D> >
   {};


} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations for base class template
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      extern template 
      class Domain<1, FFT<1>, WaveList<1>, Rpg::FieldIo<1> >;
      extern template 
      class Domain<2, FFT<2>, WaveList<2>, Rpg::FieldIo<2> >;
      extern template 
      class Domain<3, FFT<3>, WaveList<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      extern template class Domain<1>;
      extern template class Domain<2>;
      extern template class Domain<3>;
   }
}
#endif
