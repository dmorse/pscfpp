#ifndef RPC_MAX_ORDER_PARAMETER_H
#define RPC_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/MaxOrderParameterBase.h>   // base class template
#include <rpc/system/Types.h>                        // base class argument
#include <rpc/fts/analyzer/AverageAnalyzer.h>        // indirect base
#include <prdc/field/cpu/RField.h>                   // direct base member
#include <prdc/field/cpu/RFieldDft.h>                // direct base member

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Evaluates maximum squared Fourier amplitude for W_{-} field
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base template MaxOrderParameter, 
   * and inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see MaxOrderParameter
   * \see \ref rp_MaxOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class MaxOrderParameter<D, Rpc::Types<D> >
    : public MaxOrderParameterBase<D, Rpc::Types<D> >
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(Simulator<D, Rpc::Types<D> >& simulator, 
                        System<D, Rpc::Types<D> >& system);

   protected:

      /**
      * Compute and return the max order parameter.
      */
      double compute() override;

      using Base = MaxOrderParameterBase<D, Rpc::Types<D> >;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class MaxOrderParameterBase<1, Rpc::Types<1> >;
      extern template class MaxOrderParameterBase<2, Rpc::Types<2> >;
      extern template class MaxOrderParameterBase<3, Rpc::Types<3> >;
      extern template class MaxOrderParameter<1, Rpc::Types<1> >;
      extern template class MaxOrderParameter<2, Rpc::Types<2> >;
      extern template class MaxOrderParameter<3, Rpc::Types<3> >;
   }
}
#endif
