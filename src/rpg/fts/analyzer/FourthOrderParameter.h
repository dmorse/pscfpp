#ifndef RPG_FOURTH_ORDER_PARAMETER_H
#define RPG_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/FourthOrderParameterBase.h> // base template
#include <rpg/system/Types.h>                         // base argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>         // indirect base 
#include <prdc/field/cuda/RField.h>                   // base member
#include <prdc/field/cuda/RFieldDft.h>                // base member

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * FourthOrderParameter is used to detect an order-disorder transition.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp:FourthOrderParameter, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see FourthOrderDerivative
   * \see \ref rp_FourthOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class FourthOrderParameter<D, Rpg::Types<D> >
    : public FourthOrderParameterBase< D, Rpg::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(Simulator<D, Rpg::Types<D> >& simulator, 
		           System<D, Rpg::Types<D> >& system);

   private:

      /**
      * Initialize member variable prefactor_.
      */
      void computePrefactor() override;

      using Base = FourthOrderParameterBase< D, Rpg::Types<D> >;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FourthOrderParameterBase<1, Rpg::Types<1> >;
      extern template class FourthOrderParameterBase<2, Rpg::Types<2> >;
      extern template class FourthOrderParameterBase<3, Rpg::Types<3> >;
      extern template class FourthOrderParameter<1, Rpg::Types<1> >;
      extern template class FourthOrderParameter<2, Rpg::Types<2> >;
      extern template class FourthOrderParameter<3, Rpg::Types<3> >;
   }
}
#endif
