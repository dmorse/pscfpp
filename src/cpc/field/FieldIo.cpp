/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.tpp"

namespace Pscf {
   namespace Cp {
      template class FieldIo<1, CField<1,CPT>, FFT<1,CPT> >;
      template class FieldIo<2, CField<2,CPT>, FFT<2,CPT> >;
      template class FieldIo<3, CField<3,CPT>, FFT<3,CPT> >;
   }
   namespace Cpc {
      template class FieldIo<1>;
      template class FieldIo<2>;
      template class FieldIo<3>;
   } 
}
