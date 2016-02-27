#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include "Mixture.h"                       // member
#include "Grid.h"                          // member
#include <util/misc/FileMaster.h>          // member

namespace Pscf {

   class Interaction;

namespace Fd1d
{

   class Iterator;
   using namespace Util;

   /**
   * Main class in SCFT simulation of one system.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class System : public ParamComposite
   {

   public:

      /// Monomer chemical potential field type.
      typedef Propagator::WField WField;

      /// Monomer concentration / volume fraction field type.
      typedef Propagator::CField CField;

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      /// \name Lifetime (Actions)
      //@{

      /**
      * Process command line options.
      */
      void setOptions(int argc, char **argv);

      /**
      * Read input parameters (with opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      */
      void readParam();

      /**
      * Read input parameters (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read command script.
      * 
      * \param in command script file.
      */
      void readCommands(std::istream& in);

      /**
      * Read commands from default command file.
      */
      void readCommands();

      //@}
      /// \name Fields
      //@{

      /**
      * Get array of all chemical potential fields.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<WField>& wFields();

      /**
      * Get chemical potential field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      WField& wField(int monomerId);

      /**
      * Get array of all chemical potential fields.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<CField>& cFields();

      /**
      * Get chemical potential field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      CField& cField(int monomerId);

      /**
      * Read an omega field from file.
      */
      void readOmega(std::istream& in);

      /**
      * Write the current omega field to file.
      */
      void writeOmega(std::ostream& out);

      //@}
      /// \name Accessors (get objects by reference)
      //@{

      /**
      * Get Mixture by reference.
      */
      Mixture& mixture();

      /**
      * Get spatial grid by reference.
      */
      Grid& grid();

      /**
      * Get interaction by reference.
      */
      Interaction& interaction();

      /**
      * Get Iterator by reference.
      */
      Iterator& iterator();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      //@}

   private:

      // Mixture object (solves MDE for all species).
      Mixture mixture_;

      // Spatial discretization.
      Grid grid_;

      // Filemaster (holds paths to associated I/O files)
      FileMaster fileMaster_;

      // Interaction (excess free energy model)
      Interaction* interactionPtr_;

      // Pointer to associated iterator.
      Iterator* iteratorPtr_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<WField> wFields_;

      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<CField> cFields_;

      bool hasMixture_;
      bool hasGrid_;
      bool hasFields_;

      void allocateFields();

   };

   // Inline member functions

   inline Mixture& System::mixture()
   { return mixture_; }

   inline Grid& System::grid()
   { return grid_; }

   inline FileMaster& System::fileMaster()
   {  return fileMaster_; }

   inline Interaction& System::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   inline Iterator& System::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   inline 
   DArray< System::WField >& System::wFields()
   {  return wFields_; }

   inline 
   System::WField& System::wField(int id)
   {  return wFields_[id]; }

   inline
   DArray< System::CField >& System::cFields()
   {  return cFields_; }


   inline System::CField& System::cField(int id)
   {  return cFields_[id]; }

} // namespace Fd1d
} // namespace Pscf
#endif
