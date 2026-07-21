/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"

namespace Pscf {
   namespace Cp {
      using namespace Prdc::Cpu;
      template class FieldIo<1, CField<1, CppTp<1> >, FFT<1, CppTp<1> > >;
      template class FieldIo<2, CField<2, CppTp<2> >, FFT<2, CppTp<2> > >;
      template class FieldIo<3, CField<3, CppTp<3> >, FFT<3, CppTp<3> > >;
   }
   namespace Cpc {
      template class FieldIo<1>;
      template class FieldIo<2>;
      template class FieldIo<3>;
   } 
}
