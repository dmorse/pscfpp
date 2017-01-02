#ifndef UTIL_BOUNDARY_ENSEMBLE_H
#define UTIL_BOUNDARY_ENSEMBLE_H

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
   * Statistical ensemble for changes in the periodic unit cell size.
   *
   * A boundary ensemble has a type, which can be rigid or isobaric,
   * and stores a pressure if it is isobaric.
   *
   * \ingroup Ensemble_Module
   */
   class BoundaryEnsemble : public ParamComposite
   {

   public:

      /**
      * Enumeration of the allowed types of BoundaryEnsemble.
      */
      enum Type{UNKNOWN, RIGID, ISOBARIC};

      /**
      * Constructor.
      */
      BoundaryEnsemble(Type type = UNKNOWN);

      /**
      * Destructor.
      */
      ~BoundaryEnsemble();

      /**
      * Set the pressure.
      */
      void  setPressure(double pressure);

      /**
      * Read the type and (if necessary) pressure from file.
      *
      * The type is specified in the input file by a string 
      * literal "rigid" or "isobaric".
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      ///\name Accessors
      //@{
      
      /**
      * Is this an Rigid ensemble?
      */
      bool isRigid() const;

      /**
      * Is this an Isobaric ensemble?
      */
      bool isIsobaric() const;

      /**
      * Get the target pressure.
      */
      double pressure() const;

      //@}
      
      #ifdef UTIL_MPI
      /**
      * Commit associated MPI DataType.
      */
      static void commitMpiType();
      #endif

   private:

      /**
      * Target pressure 
      */
      double pressure_;

      /**
      * Ensemble type identifier.
      */
      Type  type_;

   };

   /**
   * istream extractor for an BoundaryEnsemble::Type enum value.
   *
   * \param in   input stream
   * \param type BoundaryEnsemble::Type enum value to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, BoundaryEnsemble::Type &type);

   /**
   * ostream inserter for a BoundaryEnsemble::Type enum value.
   *
   * \param  out   output stream
   * \param  type  BoundaryEnsemble::Type enum value to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const BoundaryEnsemble::Type &type);

   // Inline methods
 
   /*
   * Get the target pressure.
   */
   inline double BoundaryEnsemble::pressure() const
   {  return pressure_; }

   /*
   * Return true iff this is an rigid ensemble.
   */
   inline bool BoundaryEnsemble::isRigid() const
   {  return (type_ == RIGID); }

   /* 
   * Return true if this is an isobaric Ensemble.
   */
   inline bool BoundaryEnsemble::isIsobaric() const
   {  return (type_ == ISOBARIC); }
 
   /**
   * Serialize a BoundaryEnsemble::Type enum value.
   *
   * \param ar      archive object
   * \param data    enum value to be serialized
   * \param version archive version id
   */
   template <class Archive>
   inline void serialize(Archive& ar, BoundaryEnsemble::Type& data, const unsigned int version)
   {  serializeEnum(ar, data, version); }

   #ifdef UTIL_MPI
   /**
   * Explicit specialization MpiTraits<BoundaryEnsemble>.
   */
   template <>
   class MpiTraits<BoundaryEnsemble>
   {  
   public:  
      static MPI::Datatype type;      ///< MPI Datatype 
      static bool hasType;            ///< Is the MPI type initialized?
   };

   /**
   * Explicit specialization MpiTraits<BoundaryEnsemble::Type>.
   */
   template <>
   class MpiTraits<BoundaryEnsemble::Type>
   {  
   public:  
      static MPI::Datatype type;      ///< MPI Datatype
      static bool hasType;            ///< Is the MPI type initialized?
   };
   #endif

}
#endif
