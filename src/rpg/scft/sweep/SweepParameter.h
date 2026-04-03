#ifndef RPG_SWEEP_PARAMETER_H
#define RPG_SWEEP_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/SweepParameter.h>
#include <rpg/system/Types.h>

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;

   /**
   * Class for storing data about an individual sweep parameter.
   *
   * Specializations of this template with D=1, 2, and 3 are derived 
   * from corresponding specializations of base class template 
   * Rp::SweepParameter, and inherit their public interface and almost 
   * all of their source code from this base class.  
   *
   * \ingroup Rpg_Scft_Sweep_Module
   */
   template <int D>
   class SweepParameter
    : public Rp::SweepParameter<D, Types<D> >
   {

   public:

      /**
      * Default constructor.
      */
      SweepParameter();

      /**
      * Constructor that stores a pointer to the parent system.
      *
      * \param system  parent system
      */
      SweepParameter(System<D>& system);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SweepParameter<1, Rpg::Types<1> >;
      extern template class SweepParameter<2, Rpg::Types<2> >;
      extern template class SweepParameter<3, Rpg::Types<3> >;
      
   }
   namespace Rpg {
      extern template class SweepParameter<1>;
      extern template class SweepParameter<2>;
      extern template class SweepParameter<3>;
   }
}
#endif
