/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BoundaryEnsemble.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

namespace Util
{

   /*
   * Constructor.
   */
   BoundaryEnsemble::BoundaryEnsemble(Type type)
    : pressure_(1.0),
      type_(type)
   {  setClassName("BoundaryEnsemble"); }

   /*
   * Destructor.
   */
   BoundaryEnsemble::~BoundaryEnsemble()
   {}

   /*
   * Set the target pressure.
   */
   void BoundaryEnsemble::setPressure(double pressure)
   {
      if (!isIsobaric()) {
	 UTIL_THROW("Must be an isobaric ensemble");
      }
      pressure_ = pressure;
   }

   /*
   * Read the type and (if necessary) pressure from file.
   */
   void BoundaryEnsemble::readParameters(std::istream& in)
   {
      read<Type>(in, "type", type_);
      if (isIsobaric()) {
         read<double>(in, "pressure", pressure_);
      }
   }

   /*
   * Load internal state from an archive.
   */
   void BoundaryEnsemble::loadParameters(Serializable::IArchive &ar)
   { 
      loadParameter<Type>(ar, "type", type_);
      if (isIsobaric()) {
         loadParameter<double>(ar, "pressure", pressure_);
      }
   }

   /*
   * Save internal state to an archive.
   */
   void BoundaryEnsemble::save(Serializable::OArchive &ar)
   { 
      ar << type_;
      if (isIsobaric()) {
         ar << pressure_;
      }
   }

   /*
   * Extract an BoundaryEnsemble::Type from an istream as a string.
   */
   std::istream& operator>>(std::istream& in, BoundaryEnsemble::Type &type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "RIGID" || buffer == "rigid") {
         type = BoundaryEnsemble::RIGID;
      } else
      if (buffer == "ISOBARIC" || buffer == "isobaric") {
         type = BoundaryEnsemble::ISOBARIC;
      } else {
         UTIL_THROW("Invalid BoundaryEnsemble::Type value input");
      }
      return in;
   }

   /*
   * Insert a BoundaryEnsemble::Type to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, const BoundaryEnsemble::Type &type)
   {
      if (type == BoundaryEnsemble::RIGID) {
         out << "rigid";
      } else
      if (type == BoundaryEnsemble::ISOBARIC) {
         out << "isobaric";
      } else
      if (type == BoundaryEnsemble::UNKNOWN) {
         out << "unknown";
      }
      return out;
   }

   #ifdef UTIL_MPI
   // Initialize BoundaryEnsemble MPI Datatype.
   MPI::Datatype MpiTraits<BoundaryEnsemble>::type = MPI::BYTE;
   bool MpiTraits<BoundaryEnsemble>::hasType = false;

   // Initialize BoundaryEnsemble::Type MPI Datatype.
   MPI::Datatype MpiTraits<BoundaryEnsemble::Type>::type = MPI::INT;
   bool MpiTraits<BoundaryEnsemble::Type>::hasType = true;

   /**
   * Commit MPI Datatype.
   */
   void BoundaryEnsemble::commitMpiType()
   {
      MpiStructBuilder builder;
      BoundaryEnsemble   object;

      builder.setBase(&object);
      builder.addMember(&object.pressure_, MPI::DOUBLE);
      builder.addMember(&object.type_, MPI::INT);
      builder.commit(MpiTraits<BoundaryEnsemble>::type);
      MpiTraits<BoundaryEnsemble>::hasType = true;
   }
   #endif

}
