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

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationWriter<1, Rpg::Types<1> >;
      extern template class ConcentrationWriter<2, Rpg::Types<2> >;
      extern template class ConcentrationWriter<3, Rpg::Types<3> >;
   }
}
#endif
