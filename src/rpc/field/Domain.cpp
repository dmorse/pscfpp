/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.tpp"

namespace Pscf {
   namespace Prdc {

      using namespace Cpu;

      template class DomainReal<1, FFT<1>, WaveList<1>, Rpc::FieldIo<1> >;
      template class DomainReal<2, FFT<2>, WaveList<2>, Rpc::FieldIo<2> >;
      template class DomainReal<3, FFT<3>, WaveList<3>, Rpc::FieldIo<3> >;
   } 
   namespace Rpc {
      template class Domain<1>;
      template class Domain<2>;
      template class Domain<3>;
   } 
} 
