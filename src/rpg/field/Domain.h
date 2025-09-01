#ifndef RPG_DOMAIN_H
#define RPG_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/DomainTmpl.h>        // base class template
#include <rpg/field/FieldIo.h>            // member
#include <prdc/cuda/WaveList.h>           // member
#include <prdc/cuda/FFT.h>                // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Spatial domain for a periodic structure with real fields, on a GPU.
   *
   * See the interface of the Prdc::DomainTmpl base class template for
   * complete API documentation. The Rpg::Domain class template is 
   * basically a named partial specialization of the base class template, 
   * defined using template type parameters FFT = Prdc::Cuda::FFT<D>, 
   * WLT = Prdc::Cuda::WaveList<D>, and FIT = Rpg::FieldIo<D> . The 
   * public interface is identical to that of the base class.
   *
   * \ingroup Rpg_Field_Module
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

   // Explicit instantiation declarations
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;

} // namespace Rpg

namespace Prdc {
   // Explicit instantiation declarations for base class template
   using namespace Cuda;
   extern template 
   class DomainTmpl<1, FFT<1>, WaveList<1>, Rpg::FieldIo<1> >;
   extern template 
   class DomainTmpl<2, FFT<2>, WaveList<2>, Rpg::FieldIo<2> >;
   extern template 
   class DomainTmpl<3, FFT<3>, WaveList<3>, Rpg::FieldIo<3> >;
} // namespace Prdc

} // namespace Pscf
#endif
