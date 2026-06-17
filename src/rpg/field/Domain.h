#ifndef RPG_DOMAIN_H
#define RPG_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/Domain.h>     // base class template
#include <rpg/system/Types.h>    // base class template parameter

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
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Domain, and
   * inherit their public interface and all of their source code from
   * this base class.  
   *
   * \see Rp::Domain
   * \see \ref user_param_domain_page "Manual Page"
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Domain : public Rp::Domain< D, Types<D> >
   {
   public:
      Domain() = default;
      virtual ~Domain() = default;
   };


} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Domain<1, Rpg::Types<1> >;
      extern template class Domain<2, Rpg::Types<2> >;
      extern template class Domain<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Domain<1>;
      extern template class Domain<2>;
      extern template class Domain<3>;
   }
}
#endif
