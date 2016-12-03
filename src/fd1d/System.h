#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include "Mixture.h"                       // member
#include "Domain.h"                        // member
#include <pscf/homogeneous/Mixture.h>      // member
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter

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

      /// Base class for WField and CField
      typedef DArray<double> Field;

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

      /**
      * Compute free energy density and pressure for current fields.
      */
      void computeFreeEnergy();

      /**
      * Compute properties of homogeneous reference system.
      *
      * \int mode mode index
      */
      void computeHomogeneous(int mode);

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
      * Read chemical potential fields from file.
      *
      * \param in input stream (i.e., input file)
      */
      void readWFields(std::istream& in);

      /**
      * Write concentration or chemical potential fields to file.
      *
      * \param out output stream (i.e., output file)
      * \param fields array of fields for different species
      */
      void writeFields(std::ostream& out, Array<Field> const & fields);

      //@}
      /// \name Homogeneous reference system
      //@{

      //@}
      /// \name Accessors (get objects by reference)
      //@{

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
      * Get homogeneous mixture (used for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get the Iterator by reference.
      */
      Iterator& iterator();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      */
      double fHelmholtz() const;

      /**
      * Get precomputed pressure x monomer volume kT.
      */
      double pressure() const;

      //@}

   private:

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture mixture_;

      /**
      * Spatial domain and grid definition.
      */
      Domain domain_;

      /**
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

      /**
      * Homogeneous mixture, for reference.
      */
      Homogeneous::Mixture homogeneous_;

      /**
      * Pointer to Interaction (excess free energy model).
      */
      Interaction* interactionPtr_;

      /**
      * Pointer to associated iterator.
      */
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

      /**
      * Work array (size = # of grid points).
      */
      DArray<double> f_;

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> c_;

      /**
      * Work array (size = # of molecular species).
      */
      DArray<double> p_;

      /**
      * Work array (size = # of molecular species).
      */
      DArray<double> m_;

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Pressure times monomer volume / kT.
      */
      double pressure_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      /**
      * Have the domain and grid been initialized?
      */
      bool hasDomain_;

      /**
      * Have initial chemical potential fields been read from file?
      */
      bool hasFields_;

      void allocateFields();

      void initHomogeneous();

   };

   // Inline member functions

   /*
   * Get the associated Mixture object.
   */
   inline Mixture& System::mixture()
   { return mixture_; }

   /*
   * Get the spatial Domain.
   */
   inline Domain& System::domain()
   { return domain_; }

   /*
   * Get the FileMaster.
   */
   inline FileMaster& System::fileMaster()
   {  return fileMaster_; }

   /*
   * Get the Homogeneous::Mixture object.
   */
   inline 
   Homogeneous::Mixture& System::homogeneous()
   {  return homogeneous_; }

   /*
   * Get the Interaction (excess free energy model).
   */
   inline Interaction& System::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   /*
   * Get the Iterator (excess free energy model).
   */
   inline Iterator& System::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   /*
   * Get an array of all monomer excess chemical potential fields.
   */
   inline 
   DArray< System::WField >& System::wFields()
   {  return wFields_; }

   /*
   * Get a single monomer excess chemical potential field.
   */
   inline 
   System::WField& System::wField(int id)
   {  return wFields_[id]; }

   /*
   * Get array of all monomer concentration fields.
   */
   inline
   DArray< System::CField >& System::cFields()
   {  return cFields_; }

   /*
   * Get a single monomer concentration field.
   */
   inline System::CField& System::cField(int id)
   {  return cFields_[id]; }

   /*
   * Get precomputed Helmoltz free energy per monomer / kT.
   */
   inline double System::fHelmholtz() const
   {  return fHelmholtz_; }

   /*
   * Get precomputed pressure x monomer volume kT.
   */
   inline double System::pressure() const
   {  return pressure_; }

} // namespace Fd1d
} // namespace Pscf
#endif
