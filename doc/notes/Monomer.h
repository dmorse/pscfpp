
/*
* PSCF++ - Polymer Self-Consistent Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

   /**
   *
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
      double kuhn() const;

      /**
      * Monomer name string.
      */
      std::string name() const;

   private:

      unsigned int id_;
      double       kuhn_;
      std::string  name_;

      Field* omegaPtr_;
      Field* rhoPtr_;
   };

