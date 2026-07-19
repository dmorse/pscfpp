#ifndef RPC_FILM_ENVIRONMENT_H
#define RPC_FILM_ENVIRONMENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/*
* NOTE: In this file, we break from PSCF's conventional coding standards.
* Specifically, any subclass of MixAndMatchEnv in the Rpc namespace should 
* be declared in this header file, and the methods should be defined in 
* the class declaration rather than in a .tpp file. This is simply for the 
* sake of conciseness. These subclasses of MixAndMatchEnv require very
* little code to declare and define, so we opt to consolidate the code 
* into a single file rather than spreading it across many small files.
*/

#include "FilmFieldGenMask.h"
#include "FilmFieldGenExt.h"
#include <prdc/environment/MixAndMatchEnv.h>

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Class defining mask & external fields for thin-film systems.
   * 
   * \ingroup Rpc_Environment_Module
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
      FilmEnvironment(Rp::System<D, Rpc::Types<D> >& sys)
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
      Rp::System<D, Rpc::Types<D> >* sysPtr_;

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
         MixAndMatchEnv::fieldGenPtr1_ = new Rp::FilmFieldGenMask<D, Rpc::Types<D> >(*sysPtr_);
         MixAndMatchEnv::fieldGenPtr2_ = new Rp::FilmFieldGenExt<D, Rpc::Types<D> >(*sysPtr_);
      }

   };

   // Explicit instantiation declarations
   extern template class FilmEnvironment<1>;
   extern template class FilmEnvironment<2>;
   extern template class FilmEnvironment<3>;

} // namespace Rpc
} // namespace Pscf
#endif
