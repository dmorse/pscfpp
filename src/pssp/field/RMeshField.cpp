/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RMeshField.h"

namespace Pssp
{

   using namespace Util;

   /**
   * Default constructor.
   */
   RMeshField::RMeshField()
    : Field<double>()
   {}

   /*
   * Copy constructor.
   */
   RMeshField::RMeshField(const RMeshField& other)
    : Field<double>(other)
   {
      spaceDimension_ = other.spaceDimension_;
      meshDimensions_ = other.meshDimensions_;
   }

   /*
   * Destructor.
   */
   RMeshField::~RMeshField()
   {}

}
