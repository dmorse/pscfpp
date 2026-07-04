#ifndef RPG_MAX_ORDER_PARAMETER_H
#define RPG_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/MaxOrderParameterBase.h>   // base class template
#include <rpg/system/Types.h>                        // base class argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>        // indirect base
#include <prdc/field/cuda/RField.h>                  // direct base member
#include <prdc/field/cuda/RFieldDft.h>               // direct base member
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
   class MaxOrderParameter< D, Rpg::Types<D> >
    : public MaxOrderParameterBase<D, Rpg::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent simulator object
      * \param system  parent system object
      */
      MaxOrderParameter(Simulator<D, Rpg::Types<D> >& simulator, 
		        System<D, Rpg::Types<D> >& system);

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
      using Base = MaxOrderParameterBase<D, Rpg::Types<D> >;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class MaxOrderParameterBase<1, Rpg::Types<1> >;
      extern template class MaxOrderParameterBase<2, Rpg::Types<2> >;
      extern template class MaxOrderParameterBase<3, Rpg::Types<3> >;
      extern template class MaxOrderParameter<1, Rpg::Types<1> >;
      extern template class MaxOrderParameter<2, Rpg::Types<2> >;
      extern template class MaxOrderParameter<3, Rpg::Types<3> >;
   }
}
#endif
