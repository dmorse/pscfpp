#ifndef PRDC_SYSTEM_CONST_REF_TMPL_H
#define PRDC_SYSTEM_CONST_REF_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declarations
namespace Util {
   class FileMaster;
}

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base class template for const access to an associated System.
   *
   * A SystemConstRefTmpl holds read-only (const) pointers to an associated 
   * SystemT object (type parameter ST) and several of its primary 
   * components, for a model with real fields and periodic boundary 
   * conditions. Accessor functions return the system and its components 
   * as const references.
   *
   * \ingroup Prdc_System_Module
   */
   template <class ST>
   class SystemConstRefTmpl 
   {

   public:

      // Public type name aliases
      using SystemT = ST;
      using MixtureT = typename SystemT::MixtureT;
      using InteractionT = typename SystemT::InteractionT;
      using DomainT = typename SystemT::DomainT;
      using CFieldContainerT = typename SystemT::CFieldContainerT;
      using WFieldsT = typename SystemT::WFieldsT;
      using MaskT = typename SystemT::MaskT;
      using RFieldT = typename SystemT::RFieldT;

      // Public member functions

      /**
      * Default constructor.
      */
      SystemConstRefTmpl();

      /**
      * Constructor.
      */
      SystemConstRefTmpl(SystemT const & system);

      /**
      * Destructor.
      */
      ~SystemConstRefTmpl();

      /**
      * Create associations with a system and its components.
      */
      void associate(SystemT const & system);

      /// Get the associated System.    
      SystemT const & system() const
      {  return *systemPtr_; }

      /// Get the Mixture.
      MixtureT const & mixture() const
      {  return *mixturePtr_; }

      /// Get the Interaction.    
      InteractionT const & interaction() const
      {  return *interactionPtr_; }

      /// Get the Domain.
      DomainT const & domain() const
      {  return *domainPtr_; }

      /// Get the concentration (c) field container.
      CFieldContainerT const & c() const
      {  return *cPtr_; }

      /// Get the chemical potential (w) field container.
      WFieldsT const & w() const
      {  return *wPtr_; }

      /// Get the external potential (h) field container.
      WFieldsT const & h() const
      {  return *hPtr_; }

      /// Get the mask.
      MaskT const & mask() const
      {  return *maskPtr_; }

      /// Get the FileMaster.
      FileMaster const & fileMaster() const
      {  return *fileMasterPtr_; }

   private:

      /// Pointer to System.
      SystemT const * systemPtr_;

      /// Pointer to Mixture.
      MixtureT const * mixturePtr_;

      /// Pointer to Interaction.
      InteractionT const * interactionPtr_;

      /// Pointer to Domain.
      DomainT const * domainPtr_;

      /// Pointer to c field container.
      CFieldContainerT const * cPtr_;

      /// Pointer to chemical potential (w) field container.
      WFieldsT const * wPtr_;

      /// Pointer to external potential (h)field container.
      WFieldsT const * hPtr_;

      /// Pointer to Mask .
      MaskT const * maskPtr_;

      /// Pointer to FileMaster .
      FileMaster const * fileMasterPtr_;

   };

} // namespace Prdc
} // namespace Pscf
#endif
