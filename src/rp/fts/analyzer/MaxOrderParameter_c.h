#ifndef RPC_MAX_ORDER_PARAMETER_H
#define RPC_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/MaxOrderParameterBase.h>  // base class template
#include <pscf/backends/CPT.h>                       // base class argument

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
   class MaxOrderParameter<D,CPT>
    : public MaxOrderParameterBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(Simulator<D,CPT>& simulator, 
                        System<D,CPT>& system);

   protected:

      /**
      * Compute and return the max order parameter.
      */
      double compute() override;

      using Base = MaxOrderParameterBase<D,CPT>;

   };

   PSCF_TMPL_DECLARE_CPP(MaxOrderParameter);

} // namespace Rp
} // namespace Pscf

#if 0
// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      //extern template class MaxOrderParameterBase<1,CPT>;
      //extern template class MaxOrderParameterBase<2,CPT>;
      //extern template class MaxOrderParameterBase<3,CPT>;
      extern template class MaxOrderParameter<1,CPT>;
      extern template class MaxOrderParameter<2,CPT>;
      extern template class MaxOrderParameter<3,CPT>;
   }
}
#endif

#endif
