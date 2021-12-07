#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <fd1d/domain/Domain.h>            // member
#include <fd1d/solvers/Mixture.h>          // member
#include <pscf/homogeneous/Mixture.h>      // member
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter

namespace Pscf {

   class Interaction;

namespace Fd1d
{

   class Iterator;
   class Sweep;
   class SweepFactory;
   using namespace Util;

   /**
   * Main class in SCFT simulation of one system.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class System : public ParamComposite
   {

   public:

      /// Generic Field type.
      typedef Propagator::Field Field;

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
      *
      * This function should be called after a successful call of
      * iterator().solve(). Resulting values are returned by the 
      * freeEnergy() and pressure() accessor functions.
      */
      void computeFreeEnergy();

      /**
      * Output thermodynamic properties to a file. 
      *
      * This function outputs Helmholtz free energy per monomer,
      * pressure (in units of kT per monomer volume), and the
      * volume fraction and chemical potential of each species.
      *
      * \param out output stream 
      */
      void outputThermo(std::ostream& out);

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
      * Get the Iterator by reference.
      */
      Iterator& iterator();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double fHelmholtz() const;

      /**
      * Get precomputed pressure x monomer volume kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
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
      * Pointer to associated Sweep object
      */
      Sweep* sweepPtr_;

      /**
      * Pointer to associated Sweep object
      */
      SweepFactory* sweepFactoryPtr_;

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
      * Has the Domain been initialized?
      */
      bool hasDomain_;

      /**
      * Have initial chemical potential fields been read from file?
      */
      bool hasFields_;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep_;

      /**
      * Allocate memory for fields (private)
      */
      void allocateFields();

      /**
      * Initialize Homogeneous::Mixture object.
      */
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
   * Get precomputed pressure (units of kT / monomer volume).
   */
   inline double System::pressure() const
   {  return pressure_; }

} // namespace Fd1d
} // namespace Pscf
#endif
