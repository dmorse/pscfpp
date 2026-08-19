#ifndef RP_FILM_ENVIRONMENT_TPP
#define RP_FILM_ENVIRONMENT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/environment/FilmEnvironment.h>

namespace Pscf {
namespace Rp {

   /**
   * Create FieldGenerator objects for the mask & external field
   * 
   * This method dynamically allocates FieldGenerator objects and
   * assigns fieldGenPtr1_ and fieldGenPtr2_ to them, where the 
   * actual type of each of these objects will be a subclass of 
   * FieldGenerator, and the type will depend on the type_ parameter
   * that is read by this object.
   */
   template <int D, class T>
   void FilmEnvironment<D,T>::createGenerators()
   {
      MixAndMatchEnv::fieldGenPtr1_ = new FilmFieldGenMask<D,T>(*sysPtr_);
      MixAndMatchEnv::fieldGenPtr2_ = new FilmFieldGenExt<D,T>(*sysPtr_);
   }

} // namespace Rp
} // namespace Pscf
#endif
