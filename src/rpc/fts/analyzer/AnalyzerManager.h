#ifndef RPC_ANALYZER_MANAGER_H
#define RPC_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // direct base class template
#include <rpc/system/Types.h>                // template argument
#include <rpc/fts/analyzer/Analyzer.h>       // indirect base class member
#include <util/param/Manager.h>              // indirect base class template

// Forward declarations
namespace Pscf {
   namespace Rpc {
      template <int D> class Simulator;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Manager for a list of Analyzer objects.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::AnalyzerManager, and
   * inherit their public interface and almost all of their source code
   * from this base class.
   *
   * \see \ref rp_AnalyzerManager_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AnalyzerManager : public Rp::AnalyzerManager< D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator
      * \param system  parent System
      */
      AnalyzerManager(Simulator<D>& simulator, Rp::System<D, Rpc::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1, Rpc::Types<1> >;
      extern template class AnalyzerManager<2, Rpc::Types<2> >;
      extern template class AnalyzerManager<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class AnalyzerManager<1>;
      extern template class AnalyzerManager<2>;
      extern template class AnalyzerManager<3>;
   }
}
#endif
