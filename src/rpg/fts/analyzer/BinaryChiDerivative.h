#ifndef RPG_CHI_DERIVATIVE_H
#define RPG_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryChiDerivative.h>     // direct base template
#include <rpg/system/Types.h>                  // direct base argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>  // indirect base class

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of the base template Rp::BinaryChiDerivative, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::BinaryChiDerivative
   * \see \ref rp_BinaryChiDerivative_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryChiDerivative : public Rp::BinaryChiDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryChiDerivative(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryChiDerivative<1, Rpg::Types<1> >;
      extern template class BinaryChiDerivative<2, Rpg::Types<2> >;
      extern template class BinaryChiDerivative<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BinaryChiDerivative<1>;
      extern template class BinaryChiDerivative<2>;
      extern template class BinaryChiDerivative<3>;
   }
}
#endif
