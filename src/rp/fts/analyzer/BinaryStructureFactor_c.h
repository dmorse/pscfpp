#ifndef RPC_BINARY_STRUCTURE_FACTOR_H
#define RPC_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactorBase.h> // base template
#include <pscf/backend/cpp/CPT.h>                      // base argument
#include <pscf/backend/cpp/ConstHostArray.h>           // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Spherically averaged structure factor for a two-monomer system.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template BinaryStructureFactorBase, and inherit most
   * of their source code from this base class.
   *
   * \see BinaryStructureFactorBase
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryStructureFactor<D,CPT>
    : public BinaryStructureFactorBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryStructureFactor(Simulator<D,CPT>& simulator, 
                            System<D,CPT>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryStructureFactor<1,CPT>;
      extern template class BinaryStructureFactor<2,CPT>;
      extern template class BinaryStructureFactor<3,CPT>;
   }
}
#endif
