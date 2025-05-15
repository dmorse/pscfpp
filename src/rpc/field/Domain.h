#ifndef RPC_DOMAIN_H
#define RPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/DomainReal.h>        // base class template
#include <rpc/field/FieldIo.h>            // member
#include <prdc/cpu/WaveList.h>            // member
#include <prdc/cpu/FFT.h>                 // member
#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   template <int D>
   class Domain : 
           public Prdc::DomainReal<D, FFT<D>, WaveList<D>, FieldIo<D> >
   {

   public:

      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      ~Domain();

   };

   #ifndef RPC_DOMAIN_TPP
   // Suppress implicit instantiation
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;
   #endif

} // namespace Rpc

#ifndef RPC_DOMAIN_TPP
namespace Prdc {

   using namespace Cpu;

   extern 
   template class DomainReal<1, FFT<1>, WaveList<1>, Rpc::FieldIo<1> >;

   extern 
   template class DomainReal<2, FFT<2>, WaveList<2>, Rpc::FieldIo<2> >;

   extern 
   template class DomainReal<3, FFT<3>, WaveList<3>, Rpc::FieldIo<3> >;
} // namespace Prdc
#endif

} // namespace Pscf
#endif
