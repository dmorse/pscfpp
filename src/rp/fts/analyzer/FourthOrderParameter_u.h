#ifndef RPG_FOURTH_ORDER_PARAMETER_H
#define RPG_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/FourthOrderParameterBase.h> // base template
#include <pscf/backend/cuda/CUT.h>                        // base argument

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
   class FourthOrderParameter<D,CUT>
    : public FourthOrderParameterBase<D,CUT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(Simulator<D,CUT>& simulator, 
		           System<D,CUT>& system);

   private:

      /**
      * Initialize member variable prefactor_.
      */
      void computePrefactor() override;

      using Base = FourthOrderParameterBase<D,CUT>;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FourthOrderParameter<1,CUT>;
      extern template class FourthOrderParameter<2,CUT>;
      extern template class FourthOrderParameter<3,CUT>;
   }
}
#endif
