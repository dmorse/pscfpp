#ifndef RPC_FOURTH_ORDER_PARAMETER_H
#define RPC_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/FourthOrderParameterBase.h> // base template
#include <pscf/backends/CPT.h>                        // base argument

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * FourthOrderParameter is used to detect an order-disorder transition.
   *
   * This class evaluates the sum of fourth power of the Fourier mode 
   * amplitudes of a fluctuating exchange w field.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template FourthOrderParameterBase, and inherit their 
   * entire public interface and almost all of their source code from this
   * base class.
   *
   * \see FourthOrderParameterBase
   * \see \ref rp_FourthOrderParameter_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class FourthOrderParameter<D,CPT>
    : public FourthOrderParameterBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(
                 Simulator<D,CPT>& simulator, 
                 System<D,CPT>& system);

   private:

      /**
      * Initialize member variable prefactor_.
      */
      void computePrefactor() override;

      //  Private aliase for base class
      using Base = FourthOrderParameterBase<D,CPT>;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FourthOrderParameter<1,CPT>;
      extern template class FourthOrderParameter<2,CPT>;
      extern template class FourthOrderParameter<3,CPT>;
   }
}
#endif
