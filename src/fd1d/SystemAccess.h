#ifndef FD1D_SYSTEM_ACCESS_H
#define FD1D_SYSTEM_ACCESS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
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

      /// Base class for WField and CField
      typedef DArray<double> Field;

      /// Monomer chemical potential field type.
      typedef Propagator::WField WField;

      /// Monomer concentration / volume fraction field type.
      typedef Propagator::CField CField;

      /**
      * Constructor.
      */
      SystemAccess(System& system);

      /**
      * Destructor.
      */
      ~SystemAccess();

      /// \name Accessors (get objects by reference)
      //@{

      /**
      * Get parent System by reference.
      */
      System& system();

      /**
      * Get Mixture by reference.
      */
      Mixture& mixture();

      /**
      * Get spatial domain (including grid info) by reference.
      */
      Domain& domain();

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      Interaction& interaction();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      //@}
      /// \name Fields
      //@{

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

      //@}

   private:

      /**
      * Mixture object (solves MDE for all species).
      */ 
      System& system_;

   };

   // Inline member functions

   /*
   * Get the parent System object.
   */
   inline System& SystemAccess::system()
   { return system_; }

   /*
   * Get the associated Mixture object.
   */
   inline Mixture& SystemAccess::mixture()
   { return system_.mixture(); }

   /*
   * Get the spatial Domain.
   */
   inline Domain& SystemAccess::domain()
   { return system_.domain(); }

   /*
   * Get the FileMaster.
   */
   inline FileMaster& SystemAccess::fileMaster()
   {  return system_.fileMaster(); }

   /*
   * Get the Homogeneous::Mixture object.
   */
   inline 
   Homogeneous::Mixture& SystemAccess::homogeneous()
   {  return system_.homogeneous(); }

   /*
   * Get the Interaction (excess free energy model).
   */
   inline Interaction& SystemAccess::interaction()
   {  return system_.interaction(); }

   /*
   * Get an array of all monomer excess chemical potential fields.
   */
   inline 
   DArray< System::WField >& SystemAccess::wFields()
   {  return system_.wFields(); }

   /*
   * Get a single monomer excess chemical potential field.
   */
   inline 
   System::WField& SystemAccess::wField(int id)
   {  return system_.wField(id); }

   /*
   * Get array of all monomer concentration fields.
   */
   inline
   DArray< System::CField >& SystemAccess::cFields()
   {  return system_.cFields(); }

   /*
   * Get a single monomer concentration field.
   */
   inline System::CField& SystemAccess::cField(int id)
   {  return system_.cField(id); }

} // namespace Fd1d
} // namespace Pscf
#endif
