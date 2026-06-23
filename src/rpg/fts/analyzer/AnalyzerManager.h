#ifndef RPG_ANALYZER_MANAGER_H
#define RPG_ANALYZER_MANAGER_H

#include <rp/fts/analyzer/AnalyzerManager.h> // direct base class template
#include <rpg/system/Types.h>                // template argument
#include <rpg/fts/analyzer/Analyzer.h>       // indirect base class member
#include <util/param/Manager.h>              // indirect base class template

// Forward declarations
namespace Pscf {
   namespace Rpg {
      template <int D> class System;
      template <int D> class Simulator;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Manager for a list of Analyzer objects.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::AnalyzerManager, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::Class
   * \see \ref rp_AnalyzerManager_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
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
      AnalyzerManager(Simulator<D>& simulator, Rp::System<D, Rpg::Types<D> >& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class AnalyzerManager<1, Rpg::Types<1> >;
      extern template class AnalyzerManager<2, Rpg::Types<2> >;
      extern template class AnalyzerManager<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class AnalyzerManager<1>;
      extern template class AnalyzerManager<2>;
      extern template class AnalyzerManager<3>;
   }
}
#endif
