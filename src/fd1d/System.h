#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <fd1d/misc/FieldIo.h>             // member
#include <fd1d/solvers/Mixture.h>          // member
#include <fd1d/domain/Domain.h>            // member
#include <pscf/homogeneous/Mixture.h>      // member
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter

namespace Pscf {

   class Interaction;

namespace Fd1d
{

   class Iterator;
   class IteratorFactory;
   class Sweep;
   class SweepFactory;
   using namespace Util;

   /**
   * Main class in SCFT simulation of one system.
   *
   * A System has (among other components):
   *
   *    - a Mixture (a container for polymer and solvent solvers)
   *    - an Interaction (list of binary chi parameters)
   *    - a Domain (description of the unit cell and discretization)
   *    - monomer chemical potential fields 
   *    - monomer concentration fields 
   *    - An Iterator
   *
   * A system may also optionally contain a Sweep object.
   *
   * A minimal main program that uses this class to implement a program
   * might look something like this:
   * \code
   *   int main(int argc, char **argv) {
   *      Pscf::Fd1d::System system;
   *      system.setOptions(argc, argv);
   *      system.readParam();
   *      system.readCommands();
   *   }
   * \endcode
   * The actual main program is given in the file pscf_1d.cpp.
   *
   * \ref user_param_1d_page "Parameter File Format"
   * \ingroup Pscf_Fd1d_Module
   */
   class System : public ParamComposite
   {

   public:

      /// Generic Field type.
      typedef DArray<double> Field;

      /// Monomer chemical potential field type.
      typedef DArray<double> WField;

      /// Monomer concentration / volume fraction field type.
      typedef DArray<double> CField;

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      ///@}
      /// \name Lifetime (Actions)
      ///@{

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

      ///@}
      /// \name Primary SCFT Computations 
      //@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields (wFields).
      * This requires solution of the modified diffusion equation for 
      * all polymers, computation of Boltzmann weights for all solvents, 
      * computation of molecular partition functions for all species, 
      * computation of concentration fields for blocks and solvents, and 
      * computation of overall concentrations for all monomer types. 
      * This function does not compute the canonical (Helmholtz) free 
      * energy or grand-canonical free energy (i.e., pressure). 
      * Upon return, the flag hasCFields is set true.
      */
      void compute();
   
      /**
      * Iteratively solve a SCFT problem.
      * 
      * This function calls the iterator to attempt to solve the SCFT
      * problem for the current mixture and system parameters, using
      * the current chemical potential fields (wFields) as initial 
      * guesses.
      * 
      * \param isContinuation true if continuation within a sweep.
      * \return returns 0 for successful convergence, 1 for failure.
      */
      int iterate(bool isContinuation = false);
   
      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the 
      * parameter file to solve the SCF problem at a sequence of points
      * along a line in parameter space. The nature of this sequence of
      * points is determined by implementation of a subclass of Sweep
      * and the parameters passed to the sweep object in the parameter 
      * file.  The Iterator that is initialized in the parameter file 
      * is called at each state point.
      */
      void sweep();

      //@}
      /// \name Thermodynamic Properties
      ///@{

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * iterator().solve(). Resulting values are returned by the 
      * freeEnergy() and pressure() accessor functions.
      */
      void computeFreeEnergy();

      /**
      * Get precomputed Helmholtz free energy per monomer / kT.
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
      /// \name Thermodynamic Data Output
      ///@{

      /**
      * Write parameter file to an ostream, omitting any Sweep block. 
      *
      * \param out output stream 
      */
      void writeParamNoSweep(std::ostream& out) const;

      /**
      * Write thermodynamic properties to a file. 
      *
      * This function outputs Helmholtz free energy per monomer,
      * pressure (in units of kT per monomer volume), and the volume
      * fraction phi and chemical potential mu of each species.
      *
      * \param out output stream 
      */
      void writeThermo(std::ostream& out);

      ///@}
      /// \name Field Output 
      //@{
      
      /**
      * Write chemical potential fields in symmetrized basis format.
      *
      * \param filename name of output file
      */
      void writeW(std::string const & filename);
   
      /**
      * Write concentration fields in symmetrized basis format.
      *
      * \param filename name of output file
      */
      void writeC(std::string const & filename);
   
      /**
      * Write c-fields for all blocks and solvents in r-grid format.
      *
      * Writes concentrations for all blocks of all polymers and all
      * solvent species in r-grid format. Columns associated with blocks
      * appear ordered by polymer id and then by block id, followed by
      * solvent species ordered by solvent id. 
      *
      * \param filename name of output file
      */
      void writeBlockC(std::string  const & filename);
      
      /**
      * Write slice of a propagator at fixed s in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      * \param segmentId  integer integration step index
      */
      void writeQSlice(std::string const & filename,
                       int polymerId, int blockId,
                       int directionId, int segmentId)  const;

      /**
      * Write the final slice of a propagator in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQTail(std::string const & filename, int polymerId,
                      int blockId, int directionId)  const;

      /**
      * Write one propagator for one block, in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQ(std::string const & filename, int polymerId,
                  int blockId, int directionId)  const;

      /**
      * Write all propagators of all blocks, each to a separate file.
      *
      * Write all propagators for both directions for all blocks
      * of all polymers, with each propagator in a separate file.
      * The function writeQ is called internally for each propagator,
      * and is passed an automatically generated file name.  The file
      * name for each propagator is given by a string of the form
      * (basename)_(ip)_(ib)_(id), where (basename) denotes the value
      * of the std::string function parameter basename, and where
      * (ip), (ib), and (id) denote the string representations of
      * a polymer indiex ip, a block index ib, and direction index id,
      * with id = 0 or 1. For example, if basename == "out/q", then
      * the file name of the propagator for direction 1 of block 2
      * of polymer 0 would be "out/q_0_2_1".
      *
      * \param basename  common prefix for output file names
      */
      void writeQAll(std::string const & basename);

      //@}
      /// \name Field Accessors
      ///@{

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

      ///@}
      /// \name Member object accessors 
      ///@{

      /**
      * Get Mixture by reference.
      */
      Mixture& mixture();

      /**
      * Get Mixture by reference.
      */
      Mixture const & mixture() const;

      /**
      * Get interaction (i.e., excess free energy) by reference.
      */
      Interaction & interaction();

      /**
      * Get interaction (i.e., excess free energy) by const reference.
      */
      Interaction const & interaction() const;

      /**
      * Get spatial domain (including grid info) by reference.
      */
      Domain& domain();

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

      ///@}

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
      * FieldIo (field input-output operations).
      */
      FieldIo fieldIo_;

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
      * Pointer to associated Iterator factory.
      */
      IteratorFactory* iteratorFactoryPtr_;

      /**
      * Pointer to associated Sweep object
      */
      Sweep* sweepPtr_;

      /**
      * Pointer to associated Sweep factory.
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
      mutable DArray<double> f_;

      /**
      * Work array (size = # of monomer types).
      */
      mutable DArray<double> c_;

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Ideal gas contribution to fHelmholtz_. 
      * 
      * This encompasses the internal energy and entropy of 
      * non-interacting free chains in their corresponding 
      * potential fields defined by w_.
      */
      double fIdeal_;

      /**
      * Multi-chain interaction contribution to fHelmholtz_.
      */
      double fInter_;

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

      // Private member functions

      /**
      * Allocate memory for fields (private)
      */
      void allocateFields();

      /**
      * Initialize Homogeneous::Mixture object.
      */
      void initHomogeneous();

      /**
      * Read a string (e.g., a filename) and echo it to the log file.
      *
      * \param in  input stream from which to read
      * \param string  string variable to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

   };

   // Inline member functions

   /*
   * Get the associated Mixture object by reference.
   */
   inline Mixture& System::mixture()
   { return mixture_; }

   /*
   * Get the associated Mixture object by const reference.
   */
   inline Mixture const & System::mixture() const
   { return mixture_; }

   /*
   * Get the Interaction (excess free energy model).
   */
   inline Interaction & System::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   /*
   * Get the Interaction (excess free energy) by const reference.
   */
   inline Interaction const & System::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   /*
   * Get the spatial Domain.
   */
   inline Domain& System::domain()
   { return domain_; }

   /*
   * Get the Homogeneous::Mixture object.
   */
   inline 
   Homogeneous::Mixture& System::homogeneous()
   {  return homogeneous_; }

   /*
   * Get the Iterator.
   */
   inline Iterator& System::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   /*
   * Get the FileMaster.
   */
   inline FileMaster& System::fileMaster()
   {  return fileMaster_; }

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
