#ifndef RPG_CONCENTRATION_DERIVATIVE_H
#define RPG_CONCENTRATION_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/ConcentrationDerivative.h> // base class template
#include <rpg/system/Types.h>                        // base argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>        // indirect base

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Evaluate and average the derivative of H with respect to chi.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of base class template Rp::ConcentrationDerivative, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::ConcentrationDerivative
   * \see \ref rp_ConcentrationDerivative_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationDerivative
    : public Rp::ConcentrationDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationDerivative(Rp::Simulator<D, Rpg::Types<D> >& simulator, Rp::System<D, Rpg::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationDerivative<1, Rpg::Types<1> >;
      extern template class ConcentrationDerivative<2, Rpg::Types<2> >;
      extern template class ConcentrationDerivative<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class ConcentrationDerivative<1>;
      extern template class ConcentrationDerivative<2>;
      extern template class ConcentrationDerivative<3>;
   }
}
#endif
