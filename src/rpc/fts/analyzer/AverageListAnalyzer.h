#ifndef RPC_AVERAGE_LIST_ANALYZER_H
#define RPC_AVERAGE_LIST_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rp/fts/analyzer/AverageListAnalyzer.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates averages of several sampled real variables,
   * and optionally writes block averages to a data file during a
   * simulation. It is intended for use as a base class for Analyzers
   * that evaluate averages and (optionally) block averages for several
   * physical variables.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::AverageListAnalyzer, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::AverageListAnalyzer
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AverageListAnalyzer 
     : public Rp::AverageListAnalyzer< D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      AverageListAnalyzer(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageListAnalyzer<1, Rpc::Types<1> >;
      extern template class AverageListAnalyzer<2, Rpc::Types<2> >;
      extern template class AverageListAnalyzer<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class AverageListAnalyzer<1>;
      extern template class AverageListAnalyzer<2>;
      extern template class AverageListAnalyzer<3>;
   }
}
#endif
