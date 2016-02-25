#ifndef PSCF_SPECIES_H
#define PSCF_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
namespace Pscf
{ 

   /**
   * Base class for a molecular species (polymer or solvent).
   *
   * \ingroup Pscf_Base_Module
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
   void serialize(Archive& ar, Species::Ensemble& policy, const unsigned int version)
   { serializeEnum(ar, policy, version); }

}

#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
namespace Util
{

   /**
   * Explicit specialization MpiTraits<Species::Ensemble>.
   */
   template <>
   class MpiTraits<Pscf::Species::Ensemble> {  
   public:  
      static MPI::Datatype type;     ///< MPI Datatype
      static bool hasType;           ///< Is the MPI type initialized?
   };

}

#endif 

#endif 
