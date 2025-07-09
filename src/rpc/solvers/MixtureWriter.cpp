/*
* PSCF - MixtureWriter Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureWriter.h"
#include <prdc/solvers/MixtureWriterReal.tpp>

namespace Pscf {
   namespace Prdc { 
      template 
      class MixtureWriterReal<1, Rpc::Mixture<1>, Rpc::FieldIo<1> >;
      template 
      class MixtureWriterReal<2, Rpc::Mixture<2>, Rpc::FieldIo<2> >;
      template 
      class MixtureWriterReal<3, Rpc::Mixture<3>, Rpc::FieldIo<3> >;
   }
   namespace Rpc { 
      template class MixtureWriter<1>;
      template class MixtureWriter<2>;
      template class MixtureWriter<3>;
   }
}
