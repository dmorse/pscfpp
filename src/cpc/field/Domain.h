#ifndef CPC_DOMAIN_H
#define CPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/CppTp.h>     // backend class
#include <cp/field/Domain.h>    // base class template

// Forward declarations
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class FFT;
      template <int D, class T> class WaveList;
   }
   namespace Cpc {
      template <int D> class FieldIo;
   }
}

// Explicit instantiation declarations for base class
namespace Pscf {
   namespace Cp {
      extern template 
      class Domain<1, FFT<1, CppTp<1> >, WaveList<1, CppTp<1> >, Cpc::FieldIo<1> >;
      extern template 
      class Domain<2, FFT<2, CppTp<2> >, WaveList<2, CppTp<2> >, Cpc::FieldIo<2> >;
      extern template 
      class Domain<3, FFT<3, CppTp<3> >, WaveList<3, CppTp<3> >, Cpc::FieldIo<3> >;
   }
}

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Spatial domain for a periodic structure with real fields, on a CPU.
   *
   * The public interface of this class is identical to that of the 
   * Cp::Domain base class template. Please see documentation of
   * that base class for API documentation. 
   *
   * The Cpc::Domain\<D\> class template is a named partial specialization
   * of the base class template Cp::Domain<D, FFT, WLT, FIT> that 
   * is designed to use standard CPU hardware, defined using template type 
   * parameters FFT = Cpu::FFT\<D\>, WLT = Cpu::WaveList\<D\>,  and
   * FIT = Cpc::FieldIo\<D\> . 
   *
   * \ingroup Cpc_Field_Module
   */
   template <int D>
   class Domain 
     : public Cp::Domain< D, FFT<D, CppTp<D> >, WaveList<D, CppTp<D> >, FieldIo<D> >
   {

   public:

      /**
      * Constructor.
      *
      * Sets the class name used in the parameter file to "Domain".
      */
      Domain();

      /// Alias for base class
      using Base = Cp::Domain< D, FFT<D, CppTp<D> >, WaveList<D, CppTp<D> >, FieldIo<D> >;

      // Inherited pubic member functions
      using Base::setFileMaster;
      using Base::readParameters;
      using Base::readFieldHeader;
      using Base::unitCell;
      using Base::mesh;
      using Base::fft;
      using Base::waveList;
      using Base::fieldIo;
      using Base::lattice;

   };

   // Explicit instantiation declarations
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;

} // namespace Cpc
} // namespace Pscf
#endif
