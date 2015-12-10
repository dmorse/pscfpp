#ifndef CHEM_SPECIES_H
#define CHEM_SPECIES_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Chem
{ 

   class Species
   {
   public:
   
      enum Ensemble {UNKNOWN, CLOSED, OPEN};

      /**
      * Default constructor.
      */    
      Species();
   
      /**
      * Clear all computed quantities.
      */
      virtual void clear(){};
   
      /**
      * Solve modified diffusion equation and set related quantities.
      */
      virtual void compute(){};
   
      /**
      * Get overall volume fraction for this species.
      */
      double phi() const;
   
      /**
      * Get chemical potential for this species (units kT=1).
      */
      double mu() const;
   
      /**
      * Get molecular partition function for this species.
      */
      double q() const;
   
      /**
      * Get statistical ensemble for this species (open or closed).
      */
      Ensemble ensemble();
   
   protected:
   
      /**
      * Volume fraction, set by either setPhi or compute function.
      */
      double phi_;
   
      /**
      * Chemical potential, set by either setPhi or compute function.
      */
      double mu_;
   
      /**
      * Partition function, set by compute function.
      */
      double q_;
   
      /**
      * Statistical ensemble for this species (open or closed).
      */
      Ensemble ensemble_;
   
      /**
      * Set true by upon return by compute() and set false by clear().
      */
      bool isComputed_;
   
   };

   /*
   * Get statistical ensemble for this species (open or closed).
   */
   inline Species::Ensemble Species::ensemble()
   {  return ensemble_; }
   
} 
#endif 
