#ifndef RPC_MIXTURE_WRITER_H
#define RPC_MIXTURE_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/solvers/MixtureWriterReal.h>     // base class template
#include <rpc/solvers/Mixture.h>                // base class parameter
#include <rpc/field/FieldIo.h>                  // base class parameter

namespace Pscf {
namespace Rpc {

   /**
   * A MixtureWriter writes fields owned by a Mixture.
   *
   * An Rpc::MixWriter is a named partial specialization of
   * the class template Prdc::MixtureWriterReal. 
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class MixtureWriter 
     : public MixtureWriterReal<D, Mixture<D>, FieldIo<D> >
   {};

   // Suppress implicit instantiation
   extern template class MixtureWriter<1>;
   extern template class MixtureWriter<2>;
   extern template class MixtureWriter<3>;

} // namespace Rpc
namespace Prdc {
   // Suppress implicit instantiation of base class 
   extern template 
   class MixtureWriterReal<1, Rpc::Mixture<1>, Rpc::FieldIo<1> >;
   extern template 
   class MixtureWriterReal<2, Rpc::Mixture<2>, Rpc::FieldIo<2> >;
   extern template 
   class MixtureWriterReal<3, Rpc::Mixture<3>, Rpc::FieldIo<3> >;
} // namespace Prdc
} // namespace Pscf
#endif
