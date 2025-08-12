#ifndef RPG_FILM_ENVIRONMENT_H
#define RPG_FILM_ENVIRONMENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenMask.h"
#include "FilmFieldGenExt.h"
#include <rpg/scft/iterator/Iterator.h>
#include <prdc/environment/MixAndMatchEnv.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   template <int D> 
   class System;

   /**
   * Class defining mask & external fields for thin-film systems
   * 
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class FilmEnvironment : public MixAndMatchEnv
   {

   public:

      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      FilmEnvironment(System<D>& sys)
       : MixAndMatchEnv::MixAndMatchEnv(),
         sysPtr_(&sys)
      {  ParamComposite::setClassName("FilmEnvironment"); }

      /**
      * Destructor
      */
      ~FilmEnvironment()
      {}

   private:

      /**
      * Create FieldGenerator objects for the mask & external field
      * 
      * This method dynamically allocates FieldGenerator objects and
      * assigns fieldGenPtr1_ and fieldGenPtr2_ to them, where the 
      * actual type of each of these objects will be a subclass of 
      * FieldGenerator, and the type will depend on the type_ parameter
      * that is read by this object.
      */
      void createGenerators()
      {
         MixAndMatchEnv::fieldGenPtr1_ = new FilmFieldGenMask<D>(*sysPtr_);
         MixAndMatchEnv::fieldGenPtr2_ = new FilmFieldGenExt<D>(*sysPtr_);
      }

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Suppress implicit instantiation
   extern template class FilmEnvironment<1>;
   extern template class FilmEnvironment<2>;
   extern template class FilmEnvironment<3>;

} // namespace Rpg
} // namespace Pscf
#endif
