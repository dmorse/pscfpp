#ifndef FD1D_SYSTEM_ACCESS_H
#define FD1D_SYSTEM_ACCESS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"                        // member

namespace Pscf {
namespace Fd1d {

   using namespace Util;

   /**
   * Streamlined accesss to an associated System
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class SystemAccess 
   {

   public:

      /**
      * Default constructor.
      */
      SystemAccess();

      /**
      * Constructor.
      */
      SystemAccess(System& system);

      /**
      * Destructor.
      */
      ~SystemAccess();

      /**
      * Set the system after construction.
      */
      virtual void setSystem(System& system);

      /// \name Accessors (get objects by reference)
      ///@{

      /**
      * Get parent System by reference.
      */
      const System& system() const;

      /**
      * Get parent System by reference.
      */
      System& system();

      /**
      * Get Mixture by reference.
      */
      const Mixture& mixture() const;

      /**
      * Get Mixture by reference.
      */
      Mixture& mixture();

      /**
      * Get spatial domain (including grid info) by reference.
      */
      const Domain& domain() const;

      /**
      * Get spatial domain (including grid info) by reference.
      */
      Domain& domain();

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      const Interaction& interaction() const;

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      Interaction& interaction();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      const Homogeneous::Mixture& homogeneous() const;

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      ///@}
      /// \name Fields
      ///@{

      /**
      * Get array of all chemical potential fields.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<System::WField>& wFields();

      /**
      * Get chemical potential field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      System::WField& wField(int monomerId);

      /**
      * Get array of all chemical potential fields.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<System::CField>& cFields();

      /**
      * Get chemical potential field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      System::CField& cField(int monomerId);

      ///@}

   private:

      /**
      * Mixture object (solves MDE for all species).
      */ 
      System* systemPtr_;

   };

   // Inline member functions

   /*
   * Get the parent System object.
   */
   inline const System& SystemAccess::system() const
   { 
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_; 
   }

   /*
   * Get the parent System object.
   */
   inline System& SystemAccess::system()
   { 
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_; 
   }

   /*
   * Get the associated Mixture object.
   */
   inline const Mixture& SystemAccess::mixture() const
   { 
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->mixture(); 
   }

   /*
   * Get the associated Mixture object.
   */
   inline Mixture& SystemAccess::mixture()
   { 
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->mixture(); 
   }

   /*
   * Get the spatial Domain.
   */
   inline const Domain& SystemAccess::domain() const
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->domain(); 
   }

   /*
   * Get the spatial Domain.
   */
   inline Domain& SystemAccess::domain()
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->domain(); 
   }

   /*
   * Get the Interaction (excess free energy model).
   */
   inline const Interaction& SystemAccess::interaction() const
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->interaction(); 
   }

   /*
   * Get the Interaction (excess free energy model).
   */
   inline Interaction& SystemAccess::interaction()
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->interaction(); 
   }

   /*
   * Get the FileMaster.
   */
   inline FileMaster& SystemAccess::fileMaster()
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->fileMaster(); 
   }

   /*
   * Get the Homogeneous::Mixture object.
   */
   inline 
   Homogeneous::Mixture& SystemAccess::homogeneous()
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->homogeneous(); 
   }

   /*
   * Get an array of all monomer excess chemical potential fields.
   */
   inline 
   DArray< System::WField >& SystemAccess::wFields()
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->wFields(); 
   }

   /*
   * Get a single monomer excess chemical potential field.
   */
   inline 
   System::WField& SystemAccess::wField(int id)
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->wField(id); 
   }

   /*
   * Get array of all monomer concentration fields.
   */
   inline
   DArray< System::CField >& SystemAccess::cFields()
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->cFields(); 
   }

   /*
   * Get a single monomer concentration field.
   */
   inline System::CField& SystemAccess::cField(int id)
   {  
      UTIL_ASSERT(systemPtr_);
      return systemPtr_->cField(id); 
   }

} // namespace Fd1d
} // namespace Pscf
#endif
