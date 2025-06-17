#ifndef RPC_MIX_AND_MATCH_ENVS_H
#define RPC_MIX_AND_MATCH_ENVS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
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
#include <rpc/scft/iterator/Iterator.h>
#include <pscf/environment/MixAndMatchEnv.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> 
   class System;

   /**
   * Class defining mask & external fields for thin-film systems
   * 
   * \ingroup Rpc_Field_Module
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
      {  setClassName("FilmEnvironment"); }

      /**
      * Destructor
      */
      ~FilmEnvironment()
      {}

   protected:

      using MixAndMatchEnv::fieldGenPtr1_;
      using MixAndMatchEnv::fieldGenPtr2_;
      using ParamComposite::setClassName;

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
         fieldGenPtr1_ = new FilmFieldGenMask<D>(*sysPtr_);
         fieldGenPtr2_ = new FilmFieldGenExt<D>(*sysPtr_);
      }

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

   };

   // Suppress implicit instantiation
   extern template class FilmEnvironment<1>;
   extern template class FilmEnvironment<2>;
   extern template class FilmEnvironment<3>;

} // namespace Rpc
} // namespace Pscf
#endif
