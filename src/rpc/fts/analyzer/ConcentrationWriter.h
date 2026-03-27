#ifndef RPC_CONCENTRATION_WRITER_H
#define RPC_CONCENTRATION_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/ConcentrationWriter.h>  // base template
#include <rpc/system/Types.h>                     // template argument
#include <rpc/fts/analyzer/Analyzer.h>            // indirect base

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Periodically write c-field snapshots to a trajectory file.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::ConcentrationWriter, and inherit their 
   * entire public interface and almost all of their source code from this 
   * base class. 
   *
   * \see Rp::ConcentrationWriter
   * \see \ref rp_ConcentrationWriter_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
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
      extern template class ConcentrationWriter<1, Rpc::Types<1> >;
      extern template class ConcentrationWriter<2, Rpc::Types<2> >;
      extern template class ConcentrationWriter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class ConcentrationWriter<1>;
      extern template class ConcentrationWriter<2>;
      extern template class ConcentrationWriter<3>;
   }
}
#endif
