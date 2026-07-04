#ifndef RPC_FOURTH_ORDER_PARAMETER_H
#define RPC_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/FourthOrderParameterBase.h> // base template
#include <rpc/system/Types.h>                         // base argument
#include <rpc/fts/analyzer/AverageAnalyzer.h>         // indirect base 
#include <prdc/field/cpu/RField.h>                    // base member
#include <prdc/field/cpu/RFieldDft.h>                 // base member

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
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class FourthOrderParameter<D, Rpc::Types<D> >
    : public FourthOrderParameterBase< D, Rpc::Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(
                 Simulator<D, Rpc::Types<D> >& simulator, 
                 System<D, Rpc::Types<D> >& system);

   private:

      /**
      * Initialize member variable prefactor_.
      */
      void computePrefactor() override;

      //  Private aliase for base class
      using Base = FourthOrderParameterBase< D, Rpc::Types<D> >;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FourthOrderParameterBase<1, Rpc::Types<1> >;
      extern template class FourthOrderParameterBase<2, Rpc::Types<2> >;
      extern template class FourthOrderParameterBase<3, Rpc::Types<3> >;
      extern template class FourthOrderParameter<1, Rpc::Types<1> >;
      extern template class FourthOrderParameter<2, Rpc::Types<2> >;
      extern template class FourthOrderParameter<3, Rpc::Types<3> >;
   }
}
#endif
