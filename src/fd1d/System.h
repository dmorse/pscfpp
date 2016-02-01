#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include "Mixture.h"                       // member
#include "Grid.h"                          // member
#include <util/misc/FileMaster.h>          // member

namespace Pscf {
namespace Fd1d
{

   class Iterator;
   using namespace Util;

   class System : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      /**
      * Process command line options.
      */
      void setOptions(int argc, char **argv);

      /**
      * Read input parameters (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read input parameters (with opening and closing lines).
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      */
      void readParam();

      /**
      * Read command script.
      * 
      * \param in command script file.
      */
      void readCommands(std::istream& in);

      /**
      * Read commands from default command file.
      */
      void readCommands();

      /**
      * Get Mixture by reference.
      */
      Mixture& mixture();

      /**
      * Get spatial grid by reference.
      */
      Grid& grid();

      /**
      * Get Iterator by reference.
      */
      Iterator& iterator();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

   private:

      // Mixture object (solves MDE for all species).
      Mixture mixture_;

      // Spatial discretization.
      Grid grid_;

      // Filemaster (holds paths to associated I/O files)
      FileMaster fileMaster_;

      // Pointer to associated iterator.
      Iterator* iteratorPtr_;
   };

   // Inline member functions

   inline Mixture& System::mixture()
   { return mixture_; }

   inline Grid& System::grid()
   { return grid_; }

   inline FileMaster& System::fileMaster()
   {  return fileMaster_; }

   inline Iterator& System::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
