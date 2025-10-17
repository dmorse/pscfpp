#ifndef RPC_DOMAIN_H
#define RPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/DomainTmpl.h>        // base class template

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

// Explicit instantiation declarations of base class template
namespace Pscf {
   namespace Prdc {
      using namespace Cpu;
      extern template 
      class DomainTmpl<1, FFT<1>, WaveList<1>, Rpc::FieldIo<1> >;
      extern template 
      class DomainTmpl<2, FFT<2>, WaveList<2>, Rpc::FieldIo<2> >;
      extern template 
      class DomainTmpl<3, FFT<3>, WaveList<3>, Rpc::FieldIo<3> >;
   } 
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Spatial domain for a periodic structure with real fields, on a CPU.
   *
   * The public interface of this class is identical to that of the 
   * Prdc::DomainTmpl base class template. Please see documentation of
   * that base class for API documentation. 
   *
   * The Rpc::Domain\<D\> class template is a named partial specialization
   * of the base class template Prdc::DomainTmpl<D, FFT, WLT, FIT> that 
   * is designed to use standard CPU hardware, defined using template type 
   * parameters FFT = Prdc::Cpu::FFT\<D\>, WLT = Prdc::Cpu::WaveList\<D\>, 
   * and FIT = Rpc::FieldIo\<D\> . 
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class Domain 
     : public DomainTmpl< D, FFT<D>, WaveList<D>, FieldIo<D> >
   {

   public:

      /**
      * Constructor.
      *
      * Sets the class name used in the parameter file to "Domain".
      */
      Domain();

      /// Typedef for base class
      typedef DomainTmpl< D, FFT<D>, WaveList<D>, FieldIo<D> > Base;

      // Inherited pubic member functions
      using Base::setFileMaster;
      using Base::readParameters;
      using Base::readRGridFieldHeader;
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

   // Explicit instantiation declarations of all relevant cases
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;

} // namespace Rpc
} // namespace Pscf
#endif
