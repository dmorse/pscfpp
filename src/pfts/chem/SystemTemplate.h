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
   
   template <class TMonomer, class TPolymer, class TSolvent,
             class TWField, class TCField>
   class SystemTemplate
   {
   public:
  
      /**
      * Compute ideal gas properties for all species.
      */
      virtual void compute()
      {}
   
      /**
      * Get a Monomer type descriptor.
      * 
      * \param id integer monomer type index
      */ 
      TMonomer& monomer(int id);
  
      /**
      * Get a polymer object.
      * 
      * \param id integer polymer species index
      */ 
      TPolymer& polymer(int id);
  
      /**
      * Set a solvent solver object.
      * 
      * \param id integer solvent species index
      */ 
      TSolvent& solvent(int id);
   
      /**
      * Return W (chemical potential) field for a specific monomer type.
      * 
      * \param monomerId integer monomer type index
      */
      TWField& wField(int monomerId);
   
      /**
      * Return concentration field for specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      TCField& cField(int id);
   
      /**
      * Get number of monomers.
      */ 
      int nMonomer() const;
  
      /**
      * Get number of polymer species.
      */ 
      int nPolymer() const;
  
      /**
      * Get number of solvent (point particle) species.
      */ 
      int nSolvent() const;
  
   private:
   
      /**
      * Array of monomer type descriptors.
      */
      DArray<TMonomer> monomers_;
   
      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<TWField> wFields_;
   
      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<TCField> cFields_;
   
      /**
      * Array of polymer species solvers objects.
      *
      * Size = nPolymer.
      */
      DArray<TPolymer> polymers_;
   
      /**
      * Array of solvent species solvers.
      *
      * Size = nSolvent.
      */
      DArray<TSolvent> solvents_;
   
   };

} 
#endif 
