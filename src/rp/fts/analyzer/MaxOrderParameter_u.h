#ifndef RPG_MAX_ORDER_PARAMETER_H
#define RPG_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/MaxOrderParameterBase.h>   // base class template
#include <pscf/backends/CUT.h>                       // base class argument
#include <pscf/cuda/HostDArray.h>                    // member
#include <pscf/cuda/cudaTypes.h>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Evaluates maximum squared Fourier amplitude for W_{-} field
   *
   * \see MaxOrderParameter
   * \see \ref rp_MaxOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class MaxOrderParameter<D,CUT>
    : public MaxOrderParameterBase<D,CUT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent simulator object
      * \param system  parent system object
      */
      MaxOrderParameter(Simulator<D,CUT>& simulator, 
		        System<D,CUT>& system);

      /**
      * Setup before the start of simulation.
      */
      void setup() override;

   protected:

      /**
      * Compute and return the max order parameter.
      */
      double compute() override;

   private:

      HostDArray<cudaReal> psiHost_;

      /// Alias for base class.
      using Base = MaxOrderParameterBase<D,CUT>;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE_CUDA(MaxOrderParameter)

}
}
#endif
