/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

   /**
   * Abstract base class for a system.
   */
   class System
   {

      /**
      * Compute ideal gas properties for all species.
      */
      virtual void compute();

      /**
      * Return omega field for monomer index id.
      */
      virtual Field& omega(int id) = 0;

      /**
      * Return rho (volume fraction) field for monomer
      * index id.
      */
      virtual Field& rho(int id) = 0;

      int nMonomer();

      int nSpecies();

      Monomer& monomer(int id);

      const Array<Monomer> monomers() const;

      Species& Species(int id);

   private:

      DArray<Monomer>   monomers_;
      DPArray<Species>  species_;

   };
