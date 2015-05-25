#ifndef UTIL_SPECIES_ENSEMBLE_H
#define UTIL_SPECIES_ENSEMBLE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/archives/serialize.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>
#endif
#include <util/global.h>

namespace Util
{

   /**
   * An ensemble for the number of molecules of one Species.
   *
   * An SpeciesEnsemble has a type, which can be closed or grand. 
   * If the ensemble type is grand, it stores a chemical potential. 
   *
   * \ingroup Ensemble_Module
   */
   class SpeciesEnsemble : public ParamComposite
   {

   public:

      /**
      * Enumeration of the allowed types of SpeciesEnsemble.
      */
      enum Type{UNKNOWN, CLOSED, GRAND};

      /**
      * Constructor.
      */
      SpeciesEnsemble(Type type = UNKNOWN);

      /**
      * Destructor.
      */
      ~SpeciesEnsemble();

      /**
      * Set the chemical potential mu.
      */
      void  setMu(double mu);

      /**
      * Read the type and (if appropriate) mu from file.
      *
      * The type is specified in the input file by a string 
      * literal "closed" or "grand".
      */
      virtual void readParam(std::istream& in);

      ///\name Accessors
      //@{
      
      /**
      * Is this a Closed ensemble?
      */
      bool isClosed() const;

      /**
      * Is this a Grand ensemble?
      */
      bool isGrand() const;

      /**
      * Return the chemical potential mu.
      */
      double mu() const;

      //@}
      
      #ifdef UTIL_MPI
      /**
      * Commit associated MPI DataType.
      */
      static void commitMpiType();
      #endif

   private:

      /// Chemical potential [free energy per molecule]
      double mu_;

      /// Subclass name identifier.
      Type   type_;

   };

   /**
   * istream extractor for an SpeciesEnsemble::Type enum value.
   *
   * \param in   input stream
   * \param type SpeciesEnsemble::Type enum value to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, SpeciesEnsemble::Type &type);

   /**
   * ostream inserter for a SpeciesEnsemble::Type enum value.
   *
   * \param  out   output stream
   * \param  type  SpeciesEnsemble::Type enum value to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const SpeciesEnsemble::Type &type);

   // Inline methods
 
   /**
   * Return the mu.
   */
   inline double SpeciesEnsemble::mu() const
   {  return mu_; }

   // Return true if this is an Closed ensemble.
   inline bool SpeciesEnsemble::isClosed() const
   { return (type_ == CLOSED); }
 
   // Return true if this is an Grand Ensemble.
   inline bool SpeciesEnsemble::isGrand() const
   { return (type_ == GRAND); }

   /**
   * Serialize a SpeciesEnsemble::Type enum value.
   *
   * \param ar      archive object
   * \param data    enum value to be serialized
   * \param version archive version id
   */
   template <class Archive>
   inline void serialize(Archive& ar, SpeciesEnsemble::Type& data, const unsigned int version)
   {  serializeEnum(ar, data, version); }

   #ifdef UTIL_MPI
   /**
   * Explicit specialization MpiTraits<SpeciesEnsemble>.
   */
   template <>
   class MpiTraits<SpeciesEnsemble>
   {  
   public:  
      static MPI::Datatype type;    ///< MPI Datatype
      static bool hasType;          ///< Is the MPI type initialized?
   };

   /**
   * Explicit specialization MpiTraits<SpeciesEnsemble::Type>.
   */
   template <>
   class MpiTraits<SpeciesEnsemble::Type>
   {  
   public:  
      static MPI::Datatype type;    ///< MPI Datatype
      static bool hasType;          ///< Is the MPI type initialized?
   };
   #endif

}

#endif
