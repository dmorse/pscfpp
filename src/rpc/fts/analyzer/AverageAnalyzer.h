#ifndef RPC_AVERAGE_ANALYZER_H
#define RPC_AVERAGE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                         // grandparent base class
#include <rp/fts/analyzer/AverageAnalyzer.h>  // base class template
#include <rpc/system/Types.h>                 // class template argument

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of a single sampled real variables,
   * and optionally writes values or block averages to a data file during a
   * simulation.  It is intended for use as a base class for any Analyzer
   * that computes and evaluates an average for a single physical variable.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::AverageAnalyzer, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::AverageAnalyzer
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AverageAnalyzer : public Rp::AverageAnalyzer<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      AverageAnalyzer(Simulator<D>& simulator, System<D>& system);

      /// Alias for base class
      using AnalyzerT = typename Types<D>::Analyzer;
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   };

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AverageAnalyzer<1, Rpc::Types<1> >;
      extern template class AverageAnalyzer<2, Rpc::Types<2> >;
      extern template class AverageAnalyzer<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class AverageAnalyzer<1>;
      extern template class AverageAnalyzer<2>;
      extern template class AverageAnalyzer<3>;
   }
}
#endif
