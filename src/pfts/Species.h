#ifndef PFTS_SPECIES_H
#define PFTS_SPECIES_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
namespace Pfts
{ 

   class Species
   {
   public:
   
      enum Ensemble {Unknown, Closed, Open};

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
   class MpiTraits<Pfts::Species::Ensemble> {  
   public:  
      static MPI::Datatype type;     ///< MPI Datatype
      static bool hasType;           ///< Is the MPI type initialized?
   };

}

#endif 

#endif 
