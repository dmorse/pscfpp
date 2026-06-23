#ifndef RPC_FOURTH_ORDER_PARAMETER_H
#define RPC_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/FourthOrderParameter.h>  // base class template
#include <rpc/system/Types.h>                      // template argument
#include <rpc/fts/analyzer/AverageAnalyzer.h>      // indirect base 
#include <prdc/cpu/RField.h>                       // member
#include <prdc/cpu/RFieldDft.h>                    // member

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Prdc;

   /**
   * FourthOrderParameter is used to detect an order-disorder transition.
   *
   * This class evaluates the sum of fourth power of the Fourier mode 
   * amplitudes of a fluctuating exchange w field.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::FourthOrderParameter, and inherit their 
   * entire public interface and almost all of their source code from this
   * base class.
   *
   * \see Rp::FourthOrderParameter
   * \see \ref rp_FourthOrderParameter_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class FourthOrderParameter 
    : public Rp::FourthOrderParameter< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      FourthOrderParameter(Simulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system);

   private:

      /**
      * Initialize member variable prefactor_.
      */
      void computePrefactor() override;

      using Base = Rp::FourthOrderParameter< D, Types<D> >;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FourthOrderParameter<1, Rpc::Types<1> >;
      extern template class FourthOrderParameter<2, Rpc::Types<2> >;
      extern template class FourthOrderParameter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class FourthOrderParameter<1>;
      extern template class FourthOrderParameter<2>;
      extern template class FourthOrderParameter<3>;
   }
}
#endif
