#ifndef UTIL_ENERGY_ENSEMBLE_H
#define UTIL_ENERGY_ENSEMBLE_H

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
   * A statistical ensemble for energy.
   *
   * An energy ensemble has a type, which can be adiabatic or isothermal,
   * and stores a temperature and inverse temperature beta if the type is 
   * isothermal.
   *
   * \ingroup Ensemble_Module
   */
   class EnergyEnsemble : public ParamComposite
   {

   public:

      /**
      * Enumeration of the allowed types of EnergyEnsemble.
      */
      enum Type{UNKNOWN, ADIABATIC, ISOTHERMAL};

      /**
      * Constructor.
      */
      EnergyEnsemble(Type type = UNKNOWN);

      /**
      * Destructor.
      */
      ~EnergyEnsemble();

      /**
      * Set the temperature.
      */
      void  setTemperature(double temperature);

      /**
      * Read the type and (if necessary) temperature from file.
      *
      * The type is specified in the input file by a string 
      * literal "adiabatic" or "isothermal".
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
      * Is this an Adiabatic ensemble?
      */
      bool isAdiabatic() const;

      /**
      * Is this an Isothermal ensemble?
      */
      bool isIsothermal() const;

      /**
      * Return the temperature.
      */
      double temperature() const;

      /**
      * Return the inverse temperature.
      */
      double beta() const;

      //@}
      
      #ifdef UTIL_MPI
      /**
      * Commit MPI data type for an EnergyEnsemble.
      */
      static void commitMpiType();
      #endif

   private:

      /// Temperature * kB (units of energy)
      double temperature_;

      /// Inverse Temperature 1/(k_B T) (units of inverse energy)
      double beta_;

      /// Subclass name identifier.
      Type   type_;

   };

   /**
   * istream extractor for an EnergyEnsemble::Type enum value.
   *
   * \param in   input stream
   * \param type EnergyEnsemble::Type enum value to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, EnergyEnsemble::Type &type);

   /**
   * ostream inserter for a EnergyEnsemble::Type enum value.
   *
   * \param  out   output stream
   * \param  type  EnergyEnsemble::Type enum value to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const EnergyEnsemble::Type &type);

   // Inline methods
 
   /*
   * Return the temperature.
   */
   inline double EnergyEnsemble::temperature() const
   {  return temperature_; }

   /*
   * Return the inverse temperature.
   */
   inline double EnergyEnsemble::beta() const
   {  return beta_; }

   /*
   * Return true if this is an adiabatic ensemble.
   */
   inline bool EnergyEnsemble::isAdiabatic() const
   { return (type_ == ADIABATIC); }

   /* 
   * Return true if this is an IsothermalEnsemble.
   */
   inline bool EnergyEnsemble::isIsothermal() const
   { return (type_ == ISOTHERMAL); }

   /**
   * Serialize a EnergyEnsemble::Type enum value.
   *
   * \param ar  archive object
   * \param data  enum value to be serialized
   * \param version  archive version id
   */
   template <class Archive>
   inline void serialize(Archive& ar, EnergyEnsemble::Type& data, const unsigned int version)
   {  serializeEnum(ar, data, version); }


   #ifdef UTIL_MPI
   /**
   * Explicit specialization MpiTraits<EnergyEnsemble>.
   */
   template <>
   class MpiTraits<EnergyEnsemble>
   {  
   public:  
      static MPI::Datatype type;       ///< MPI Datatype
      static bool hasType;             ///< Is the MPI type initialized?  
   };

   /**
   * Explicit specialization MpiTraits<EnergyEnsemble::Type>.
   */
   template <>
   class MpiTraits<EnergyEnsemble::Type>
   {  
   public:  
      static MPI::Datatype type;       ///< MPI Datatype
      static bool hasType;             ///< Is the MPI type initialized?
   };
   #endif

}
#endif
