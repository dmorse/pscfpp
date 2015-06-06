#ifndef PFTS_SYSTEM_TEMPLATE_H
#define PFTS_SYSTEM_TEMPLATE_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{ 

   using namespace Util;
   
   template <class Monomer, class Polymer, class Solvent,
             class WField, class CField>
   class SystemTemplate
   {
   public:
  
      /**
      * Get number of monomers.
      */ 
      int nMonomer();
  
      /**
      * Get number of polymer species.
      */ 
      int nPolymer();
  
      /**
      * Get number of solvent (point particle) species.
      */ 
      int nSolvent();
  
      /**
      * Get a Monomer type descriptor.
      * 
      * \param id integer monomer type index
      */ 
      Monomer& monomer(int id);
  
      /**
      * Set a Polymer solver object.
      * 
      * \param id integer polymer species index
      */ 
      Polymer& polymer(int id);
  
      /**
      * Set a Solvent solver object.
      * 
      * \param id integer solvent species index
      */ 
      Solvent& solvent(int id);
   
      /**
      * Compute ideal gas properties for all species.
      */
      virtual void compute()
      {}
   
      /**
      * Return W (chemical potential) field for monomer index id.
      */
      WField& wField(int id);
   
      /**
      * Return concentration field for monomer index id.
      */
      CField& cField(int id);
   
   private:
   
      /**
      * Array of monomer type descriptors.
      */
      DArray<Monomer> monomers_;
   
      /**
      * Array of fields associated with monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<WField> wFields_;
   
      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<CField> cFields_;
   
      /**
      * Array of polymer species solvers objects.
      */
      DArray<Polymer> polymers_;
   
      /**
      * Array of solvent species solvers.
      */
      DArray<Solvent> solvents_;
   
   };

} 
#endif 
