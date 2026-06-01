#ifndef RPG_CUBIC_LENGTH_DERIVATIVE_H
#define RPG_CUBIC_LENGTH_DERIVATIVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/CubicLengthDerivative.h> // direct base template
#include <rpg/system/Types.h>                      // direct base argument
#include <rpg/fts/analyzer/AverageAnalyzer.h>      // indirect base class

namespace Pscf {
namespace Rpg {

   // Forward declarations
   template <int D> class System;
   template <int D> class Simulator;

   /**
   * Evaluate the derivative of H with respect to cubic box length.
   *
   * This class is designed specifically for use with a cubic unit cell.
   *
   * Specializations of this template with D=1, 2 and 3 are derived from
   * specializations of the base class template Rp::CubicLengthDerivative,
   * and inherit their public interface and almost all of their source
   * code from this base class.
   *
   * \see Rp::CubicLengthDerivative
   * \see \ref rp_CubicLengthDerivative_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class CubicLengthDerivative
     : public Rp::CubicLengthDerivative< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      */
      CubicLengthDerivative(Simulator<D>& simulator, System<D>& system);

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class CubicLengthDerivative<1, Rpg::Types<1> >;
      extern template class CubicLengthDerivative<2, Rpg::Types<2> >;
      extern template class CubicLengthDerivative<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class CubicLengthDerivative<1>;
      extern template class CubicLengthDerivative<2>;
      extern template class CubicLengthDerivative<3>;
   }
}
#endif
