#ifndef PSCF_SPECIES_H
#define PSCF_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
namespace Pscf
{ 

   /**
   * Base class for a molecular species (polymer or solvent).
   *
   * \ingroup Pscf_Chem_Module
   */
   class Species
   {
   public:
  
      /**
      * Statistical ensemble for number of molecules.
      */ 
      enum Ensemble {Unknown, Closed, Open};

      /**
      * Default constructor.
      */    
      Species();
   
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
   
   };

   /*
   * Get species volume fraction.
   */
   inline double Species::phi() const
   {  return phi_; }
   
   /*
   * Get species chemical potential.
   */
   inline double Species::mu() const
   {  return mu_; }
   
   /*
   * Get statistical ensemble for this species (open or closed).
   */
   inline Species::Ensemble Species::ensemble()
   {  return ensemble_; }
   
   /**
   * istream extractor for a Species::Ensemble.
   *
   * \param  in       input stream
   * \param  policy   Species::Ensemble to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Species::Ensemble& policy);

   /**
   * ostream inserter for an Species::Ensemble.
   *
   * \param  out      output stream
   * \param  policy   Species::Ensemble to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, Species::Ensemble policy);

   /**
   * Serialize a Species::Ensemble
   *
   * \param ar      archive object
   * \param policy  object to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, Species::Ensemble& policy, 
                  const unsigned int version)
   { serializeEnum(ar, policy, version); }

}
#endif 
