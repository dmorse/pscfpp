#ifndef PSCF_CHEM_MONOMER_H
#define PSCF_CHEM_MONOMER_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <string>

namespace Pfts{ 

   /**
   * Descriptor for a monomer or particle type.
   */
   class Monomer
   {
   public:

      /**
      * Constructor.
      */
      Monomer();

      /**
      * Unique integer index for monomer type.
      */
      int id() const;

      /**
      * Statistical segment length.
      */
      double step() const;

      /**
      * Monomer name string.
      */
      std::string name() const;

   private:

      int  id_;
      double  step_;
      std::string  name_;

   };

} 
#endif 
