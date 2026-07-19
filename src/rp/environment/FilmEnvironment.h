#ifndef RP_FILM_ENVIRONMENT_H
#define RP_FILM_ENVIRONMENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/environment/MixAndMatchEnv.h>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class FilmFieldGenMask;
      template <int D, class T> class FilmFieldGenExt;
      template <int D, class T> class System;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Class defining mask & external fields for thin-film systems.
   * 
   * \ingroup Rp_Environment_Module
   */
   template <int D, class T>
   class FilmEnvironment : public MixAndMatchEnv
   {

   public:

      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      FilmEnvironment(System<D,T>& sys)
       : MixAndMatchEnv::MixAndMatchEnv(),
         sysPtr_(&sys)
      {  ParamComposite::setClassName("FilmEnvironment"); }

      /**
      * Destructor
      */
      ~FilmEnvironment()
      {}

   private:

      /// Pointer to the associated system object.
      System<D,T>* sysPtr_;

      /**
      * Create FieldGenerator objects for the mask & external field
      * 
      * This method dynamically allocates FieldGenerator objects and
      * assigns fieldGenPtr1_ and fieldGenPtr2_ to them, where the 
      * actual type of each of these objects will be a subclass of 
      * FieldGenerator, and the type will depend on the type_ parameter
      * that is read by this object.
      */
      void createGenerators();

   };

} // namespace Rp
} // namespace Pscf
#endif
