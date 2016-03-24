#ifndef PSCF_UNIT_CELL_H
#define PSCF_UNIT_CELL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCellTmpl.h"

namespace Pscf
{ 

   using namespace Util;

   /**
   * Crystal unit cell.
   *
   * \ingroup Pscf_Base_Module
   */
   template <int D>
   class UnitCell : public UnitCellTmpl<D>
   {
   public:

      /**
      * Constructor.
      */
      UnitCell(){}

   
      /**
      * Destructor.
      */
      ~UnitCell(){}
   
   };
   
} 
#endif 
