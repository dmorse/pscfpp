#ifndef PSCF_CHEM_MONOMER_H
#define PSCF_CHEM_MONOMER_H

/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf{ 
namespace Chem{

   /**
   * Descriptor for a monomer or particle type.
   */
   class Monomer
   {
   public:

      /**
      * Unique integer index for monomer type.
      */
      unsigned int id() const;

      /**
      * Statistical segment length.
      */
      double step() const;

      /**
      * Monomer name string.
      */
      std::string name() const;

   private:

      unsigned int  id_;
      double  step_;
      std::string  name_;

   };

} 
} 
#endif 
