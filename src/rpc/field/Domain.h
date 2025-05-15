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

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Spatial domain for a periodic structure with real fields, on a CPU.
   *
   * See the interface of the Prdc::DomainReal base class template for
   * complete interface documentation. The Rpc::Domain class template is 
   * basically a named partial specialization of the base class template, 
   * defined using template type parameters FFT = Prdc::Cpu::FFT<D>, 
   * WLT = Prdc::Cpu::WaveList<D>, and FIT = Rpc::FieldIo<D> . The 
   * public interface is identical to that of the base class.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Domain 
     : public DomainReal< D, FFT<D>, WaveList<D>, FieldIo<D> >
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

      /// Typedef for base class
      typedef DomainReal< D, FFT<D>, WaveList<D>, FieldIo<D> > Base;

      // Inherited pubic member functions

      using Base::setFileMaster;
      using Base::readParameters;
      using Base::readRGridFieldHeader;
      using Base::setUnitCell;
      using Base::makeBasis;
      using Base::unitCell;
      using Base::mesh;
      using Base::group;
      using Base::basis;
      using Base::fft;
      using Base::waveList;
      using Base::fieldIo;
      using Base::lattice;
      using Base::groupName;
      using Base::hasGroup;
      using Base::hasBasis;
      using Base::writeStars;
      using Base::writeWaves;
      using Base::writeGroup;

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
