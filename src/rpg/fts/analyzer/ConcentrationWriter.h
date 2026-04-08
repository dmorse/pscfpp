#ifndef RPG_CONCENTRATION_WRITER_H
#define RPG_CONCENTRATION_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/ConcentrationWriter.h>  // base template
#include <rpg/system/Types.h>                     // template argument
#include <rpg/fts/analyzer/Analyzer.h>            // indirect base

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Periodically write c-field snapshots to a trajectory file.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::ConcentrationWriter, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::ConcentrationWriter
   * \see \ref rp_ConcentrationWriter_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class ConcentrationWriter
    : public Rp::ConcentrationWriter<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      ConcentrationWriter(Simulator<D>& simulator, System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationWriter<1, Rpg::Types<1> >;
      extern template class ConcentrationWriter<2, Rpg::Types<2> >;
      extern template class ConcentrationWriter<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class ConcentrationWriter<1>;
      extern template class ConcentrationWriter<2>;
      extern template class ConcentrationWriter<3>;
   }
}
#endif
