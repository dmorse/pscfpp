#ifndef RP_SYSTEM_CONST_REF_H
#define RP_SYSTEM_CONST_REF_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declaration
namespace Util {
   class FileMaster;
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

      // Public type name aliases
      using SystemT = typename T::System;
      using MixtureT = typename T::Mixture;
      using InteractionT = typename T::Interaction;
      using DomainT = typename T::Domain;
      using CFieldsT = typename T::CFields;
      using WFieldsT = typename T::WFields;
      using MaskT = typename T::Mask;
      using RFieldT = typename T::RField;

      // Public member functions

      /**
      * Default constructor.
      */
      SystemConstRef();

      /**
      * Constructor that creates associations.
      * 
      * Using of this constructor is equivalent to using the default
      * constructor and then invoking the associate() function.
      * 
      * \param system  parent system object
      */
      SystemConstRef(typename T::System const & system);

      /**
      * Destructor.
      */
      ~SystemConstRef();

      /**
      * Create associations with a system and its components.
      * 
      * \param system  parent system object
      */
      void associate(typename T::System const & system);

      /**
      * Get the associated System.
      */
      typename T::System const & system() const
      {  return *systemPtr_; }

      /**
      * Get the Mixture.
      */
      typename T::Mixture const & mixture() const
      {  return *mixturePtr_; }

      /**
      * Get the Interaction.
      */
      typename T::Interaction const & interaction() const
      {  return *interactionPtr_; }

      /**
      * Get the Domain.
      */
      typename T::Domain const & domain() const
      {  return *domainPtr_; }

      /**
      * Get the concentration (c) field container.
      */
      typename T::CFields const & c() const
      {  return *cPtr_; }

      /**
      * Get the chemical potential (w) field container.
      */
      typename T::WFields const & w() const
      {  return *wPtr_; }

      /**
      * Get the external potential (h) field container (if any).
      */
      typename T::WFields const & h() const
      {  return *hPtr_; }
     
      /** 
      * Get the mask (if any).
      */
      typename T::Mask const & mask() const
      {  return *maskPtr_; }

      /**
      * Get the FileMaster.
      */
      FileMaster const & fileMaster() const
      {  return *fileMasterPtr_; }

   private:

      /// Pointer to System.
      typename T::System const * systemPtr_;

      /// Pointer to Mixture.
      typename T::Mixture const * mixturePtr_;

      /// Pointer to Interaction.
      typename T::Interaction const * interactionPtr_;

      /// Pointer to Domain.
      typename T::Domain const * domainPtr_;

      /// Pointer to monomer concentration (c) field container.
      typename T::CFields const * cPtr_;

      /// Pointer to chemical potential (w) field container.
      typename T::WFields const * wPtr_;

      /// Pointer to external potential (h) field container.
      typename T::WFields const * hPtr_;

      /// Pointer to Mask.
      typename T::Mask const * maskPtr_;

      /// Pointer to FileMaster .
      FileMaster const * fileMasterPtr_;

   };

} // namespace Rp
} // namespace Pscf
#endif
