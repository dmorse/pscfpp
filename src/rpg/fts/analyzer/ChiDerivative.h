#ifndef RPG_CHI_DERIVATIVE_H
#define RPG_CHI_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                // indirect base class
#include <rp/fts/analyzer/ChiDerivative.h>  // base class template
#include <rpg/system/Types.h>               // base template argument

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate the derivative of H with respect to chi.
   *
   * Specializations of this template are derived from specializations of
   * the base class template Rp::ChiDerivative, and inherit their entire
   * public interface and almost all of their source code from this base
   * class.
   *
   * \see Rp::ChiDerivative
   * \see \ref rp_ChiDerivative_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class ChiDerivative : public Rp::ChiDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ChiDerivative(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ChiDerivative<1, Rpg::Types<1> >;
      extern template class ChiDerivative<2, Rpg::Types<2> >;
      extern template class ChiDerivative<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class ChiDerivative<1>;
      extern template class ChiDerivative<2>;
      extern template class ChiDerivative<3>;
   }
}
#endif
