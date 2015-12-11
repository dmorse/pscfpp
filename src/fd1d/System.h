#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pfts/SystemTmpl.h>

namespace Pfts {
namespace Fd1d
{

   class System : public SystemTmpl<Polymer, Solvent>
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
      * Read all parameters and initialize.
      */
      void readParameters(std::istream& in);

      /**
      * Compute molecular partition functions and concentrations.
      */
      void compute();

      /**
      * Get number of spatial grid points.
      */
      int nx() const;

      /**
      * Get spatial grid step size.
      */
      double dx() const;

   private:

      double xMin_;
      double xMax_;
      double dx_;
      double ds_;
      int nx_;

   };

   // Inline member functions

   inline int System::nx() const
   {  return nx_; }

   inline double System::dx() const
   {  return dx_; }

}
}
#endif
