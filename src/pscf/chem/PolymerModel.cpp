/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerModel.h"
#include <util/global.h>
#include <string>

// Anonymous namespace containing pseudo-private variables.
namespace {

  // Comments:
  //
  //  1) Variables model_, nSet_ and isLocked_ defined in this anonmymous
  //     namespace are only accessible in this file, and are thus 
  //     pseudo-private
  //
  //  2) model_ is initialized to PolymerModel::Thread
  //  
  //  3) If isLocked_ , model_ may never be altered again.
  //
  //  4) isLocked_ is initialized to false

  Pscf::PolymerModel::Type model_ = Pscf::PolymerModel::Thread;

  bool isLocked_ = false;

  int nSet_ = 0;

}

namespace Pscf { 
namespace PolymerModel {
   
   using namespace Util;

   // Mutators

   /*
   * Set the global polymer model enumeration value.
   */
   void setModel(Type model)
   {  
      UTIL_CHECK(!isLocked_);
      model_ = model; 
      ++nSet_;
   }

   /*
   * Permanently lock the model type.
   */
   void lock() 
   {  isLocked_ = true; }

   // Accessors

   /*
   * Get the global polymer model enumeration value.
   */
   Type model()
   {  return model_; }

   /*
   * Is the global polymer model a continuous thread model?
   */
   bool isThread()
   {  return (model_ == Type::Thread); }

   /*
   * Is the global polymer model a discrete bead model?
   */
   bool isBead()
   {  return (model_ == Type::Bead); }

   /*
   * Is the model type locked?
   */
   bool isLocked() 
   {  return isLocked_; }

   /*
   * How many times has the choice of model been set ? 
   */
   int nSet() 
   {  return nSet_; }
   
   /*
   * Input stream extractor for a PolymerModel::Type enumeration.
   */ 
   std::istream& operator >> (std::istream& in, 
                              PolymerModel::Type& type)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Thread" || buffer == "thread") {
         type = PolymerModel::Thread;
      } else
      if (buffer == "Bead" || buffer == "bead") {
         type = PolymerModel::Type::Bead;
      } else {
         std::string msg = "Unknown input PolymerModel value string: ";
         msg += buffer;
         UTIL_THROW(msg.c_str());
      } 
      return in;
   }

   /*
   * Output stream inserter for a PolymerModel::Type enumeration.
   */ 
   std::ostream& operator << (std::ostream& out, 
                              PolymerModel::Type const & type)
   {
      if (type == PolymerModel::Thread) {
         out << "thread";
      } else
      if (type == PolymerModel::Type::Bead) {
         out << "bead";
      } else {
         UTIL_THROW("Error writing a PolymerModel value");
      }
      return out;
   }

} // end namespace PolymerModel
} // end namespace Pscf
