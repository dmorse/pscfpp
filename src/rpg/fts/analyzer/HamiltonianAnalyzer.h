#ifndef RPG_HAMILTONIAN_ANALYZER_H
#define RPG_HAMILTONIAN_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"                  // indirect base class
#include <rp/fts/analyzer/HamiltonianAnalyzer.h>  // base class template
#include <rpg/system/Types.h>                     // template argument

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Compute averages and output block averages of Hamiltonian components.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::HamiltonianAnalyzer , and inherit their 
   * entire public interface and almost all of their source code from this 
   * base class. 
   *
   * \see \ref rp_HamiltonianAnalyzer_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class HamiltonianAnalyzer 
     : public Rp::HamiltonianAnalyzer< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      HamiltonianAnalyzer(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class HamiltonianAnalyzer< 1, Rpg::Types<1> >;
      extern template class HamiltonianAnalyzer< 2, Rpg::Types<2> >;
      extern template class HamiltonianAnalyzer< 3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class HamiltonianAnalyzer<1>;
      extern template class HamiltonianAnalyzer<2>;
      extern template class HamiltonianAnalyzer<3>;
   }
}
#endif
