#ifndef RPG_PERTURBATION_DERIVATIVE_H
#define RPG_PERTURBATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <rp/fts/analyzer/PerturbationDerivative.h>
#include <rpg/system/Types.h>

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Evaluate derivative of H w/ respect to perturbation parameter lambda.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of base class template Rp::PerturbationDerivative, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::PerturbationDerivative
   * \see rp_PerturbationDerivative_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class PerturbationDerivative 
     : public Rp::PerturbationDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      PerturbationDerivative(Simulator<D>& simulator, Rp::System<D, Rpg::Types<D> >& system);

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PerturbationDerivative< 1, Rpg::Types<1> >;
      extern template class PerturbationDerivative< 2, Rpg::Types<2> >;
      extern template class PerturbationDerivative< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class PerturbationDerivative<1>;
      extern template class PerturbationDerivative<2>;
      extern template class PerturbationDerivative<3>;
   }
}
#endif
