#ifndef RP_SYSTEM_CONST_REF_H
#define RP_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/TmplDeclare.h>   // declaration macros

// Forward declarations
namespace Util {
   class FileMaster;
}
namespace Pscf {
   class Interaction;
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Mixture;
      template <int D, class T> class Domain;
      template <int D, class T> class CFields;
      template <int D, class T> class WFields;
      template <int D, class T> class Mask;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Base class template for const access to an associated System.
   *
   * A SystemConstRef holds read-only (const) pointers to an associated 
   * System object (type T::System) and several of its primary components,
   * for a model with real fields and periodic boundary conditions. 
   * Accessor functions return the system and its components as const 
   * references.
   *
   * \ingroup Rp_System_Module
   */
   template <int D, class T>
   class SystemConstRef 
   {

   public:

      /**
      * Get the associated System.
      */
      System<D,T> const & system() const
      {  return *systemPtr_; }

      /**
      * Get the Mixture.
      */
      Mixture<D,T> const & mixture() const
      {  return *mixturePtr_; }

      /**
      * Get the Interaction.
      */
      Interaction const & interaction() const
      {  return *interactionPtr_; }

      /**
      * Get the Domain.
      */
      Domain<D,T> const & domain() const
      {  return *domainPtr_; }

      /**
      * Get the concentration (c) field container.
      */
      CFields<D,T> const & c() const
      {  return *cPtr_; }

      /**
      * Get the chemical potential (w) field container.
      */
      WFields<D,T> const & w() const
      {  return *wPtr_; }

      /**
      * Get the external potential (h) field container (if any).
      */
      WFields<D,T> const & h() const
      {  return *hPtr_; }
     
      /** 
      * Get the mask (if any).
      */
      Mask<D,T> const & mask() const
      {  return *maskPtr_; }

      /**
      * Get the FileMaster.
      */
      FileMaster const & fileMaster() const
      {  return *fileMasterPtr_; }

   protected:

      /**
      * Default constructor.
      */
      SystemConstRef();

      /**
      * Constructor - creates associations.
      * 
      * Using of this constructor is equivalent to using the default
      * constructor and then invoking the associate() function.
      * 
      * \param system  parent system object
      */
      SystemConstRef(System<D,T> const & system);

      /**
      * Destructor.
      */
      ~SystemConstRef();

      /**
      * Create associations with a system and its components.
      * 
      * \param system  parent system object
      */
      void associate(System<D,T> const & system);

   private:

      /// Pointer to System.
      System<D,T> const * systemPtr_;

      /// Pointer to Mixture.
      Mixture<D,T> const * mixturePtr_;

      /// Pointer to Interaction.
      Interaction const * interactionPtr_;

      /// Pointer to Domain.
      Domain<D,T> const * domainPtr_;

      /// Pointer to monomer concentration (c) field container.
      CFields<D,T> const * cPtr_;

      /// Pointer to chemical potential (w) field container.
      WFields<D,T> const * wPtr_;

      /// Pointer to external potential (h) field container.
      WFields<D,T> const * hPtr_;

      /// Pointer to Mask.
      Mask<D,T> const * maskPtr_;

      /// Pointer to FileMaster .
      FileMaster const * fileMasterPtr_;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(SystemConstRef)

} // namespace Rp
} // namespace Pscf
#endif
