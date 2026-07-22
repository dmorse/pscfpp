#ifndef CP_SYSTEM_H
#define CP_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <pscf/chem/PolymerModel.h>      // member
#include <pscf/backends/CPT.h>              // backend class

// Forward declarations
namespace Util {
   class FileMaster;
   template <typename E, int N> class FSArray;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
      template <int D, class T> class CField;
   }
}

namespace Pscf {
namespace Cp {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * Base class template for a complex CL-FTS system.
   *
   * <b> Template parameters and typename aliases</b>:
   *
   *    D - integer dimensionality of space (D=1, 2, or 3)
   *    T - "Types" class collection of aliases for related classes
   * 
   * <b> Usage </b>: A specialization of System\<D, T\> is used as a base 
   * class for System\<D\> classes defined in namespaces Cpc and Cpg, 
   * for D=1, 2, or 3.  In this usage, template parameter T must be a
   * specialization named Types<D> of a template \<int D\> class Types that
   * is defined in each of these two namespaces. For example, in namespace
   * Cpc, for each value of D, class Cpc::System\<D\> is derived from the
   * class Prdc::System\< D, Cpc::Types\<D\> >. For each such program
   * level namespace, Types\<D\> defines a set of of typename aliases 
   * for classes used in that namespace, for the specified value of D. 
   * For example, the typename Cpc::Types\<D\>::Mixture is an alias for 
   * the type Cpc::Mixture<D> used to represent a mixture in the Cpc 
   * namespace for systems of dimension D. See the definitions of 
   * Cpc::Types and Cpg::Types for details. 
   *
   * In the remainder of the documentation for this template, Cp::System, 
   * unqualified names such as "Mixture", "Iterator", etc. are often used 
   * as shorthand for typename aliases such as T::Mixture, T::Iterator 
   * that are defined in class T (i.e., in Cpc::Types\<D\> or 
   * Cpg::Types\<D\>)
   *
   * <b> Class Components </b>:
   * A System object has:
   *
   *    - a Mixture (container for polymer and solvent solvers)
   *    - an %Interaction (description of binary interactions)
   *    - a Domain (description of unit cell and discretization)
   *    - a WFields container of monomer chemical potential (w) fields
   *    - a CFields container of monomer concentration (c) fields
   *    - a Simulator
   *
   * \ingroup Cp_System_Module
   */
   template <int D, class T>
   class System : public ParamComposite
   {

   public:

      // Public type name aliases
      using MixtureT = typename T::Mixture;
      using InteractionT = typename T::Interaction;
      using DomainT = typename T::Domain;
      using WFieldsT = typename T::WFields;  // chemical pot. fields
      using CFieldsT = typename T::CFields;  // concentration fields
      using CFieldT = CField<D,CPT>;    // complex field type

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      *
      * When a specialization of System<D,T> is used as a base class for 
      * a concrete System class, such as Cpc::System\<D\>, the typename 
      * T::System must be an alias for the name of the subclass. In 
      * this usage, in the member initialization list of the T::System 
      * subclass constructor,  a reference to the subclass instance must
      * be passed as "*this" to this System base class constructor. 
      * The address of the instance of the T::System subclass is then 
      * retained in the Cp::System base class instance by a private 
      * member variable named systemPtr_ that is of type T::System* . 
      * See definitions of the constructors for the Cpc::System and 
      * Cpc::System class templates for examples of this usage.
      *
      * \param system  instance of System subclass
      */
      System(typename T::System& system);

      /**
      * Destructor.
      */
      ~System();

      // Suppress several compiler-generated member functions
      System() = delete;
      System(System<D,T> const &) = delete;
      System<D,T>& operator = (System<D,T> const & ) = delete;

      ///@}
      /// \name Lifetime Actions 
      ///@{

      /**
      * Process command line options.
      *
      * This function takes the same arguments as any C/C++ main program
      * function. The arguments of the main function should d be passed
      * to this function unaltered, to allow this function to process 
      * the command line options.
      *
      * \param argc  number of command line arguments
      * \param argv  array of pointers to command line arguments
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
      *
      * This function reads the parameter file set by the -p command
      * line option.
      */
      void readParam();

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read and process commands from an input stream.
      *
      * \param in command script input stream
      */
      void readCommands(std::istream& in);

      /**
      * Read and process commands from the default command file.
      *
      * This function reads the parameter file set by the -c command
      * line option.
      */
      void readCommands();

      ///@}
      /// \name Field Theory Computations
      ///@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent system chemical potential fields. This
      * requires solution of the modified diffusion equation for all
      * polymers, computation of Boltzmann weights for all solvents,
      * computation of molecular partition functions for all species,
      * computation of concentration fields for all blocks and solvents,
      * and computation of overall concentrations for all monomer types.
      * This function does not compute the canonical (Helmholtz) free
      * energy or grand-canonical free energy (i.e., pressure).
      *
      * This function also computes the stress, by calling 
      * computeStress(), if and only if the argument needStress is true. 
      *
      * \pre  w().hasData() == true
      * \post c().hasData() == true
      * \post hasStress() == true iff needStress == true
      *
      * \param needStress  true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Perform a field theoretic simulation (PS-FTS).
      *
      * Perform a complex Langevin field theoretic simulation (Cl-FTS).
      * The number of CL steps to be performed is given by the parameter
      * "nStep".
      *
      * \pre Function w().hasData() == true
      *
      * \param nStep  number of simulation (BD or MC) steps
      */
      void simulate(int nStep);

      /**
      * Mark c-fields and free energy as outdated or invalid.
      *
      * This function should be called whenever any of the inputs to the
      * solution of the modified diffusion equation are modified, including
      * the w fields or unit cell parameters. Upon return, c().hasData() 
      * will return false.
      */
      void clearCFields();

      ///@}
      /// \name Unit Cell Modifiers
      ///@{

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice (i.e., lattice system type) set in the UnitCell<D>
      * unitCell input parameter must agree with the lattice enum value
      * that was set previously in the parameter file.
      *
      * If a space group has been declared but a basis has not yet been
      * initialized, then a symmetry-adapted basis will be constructed.
      *
      * \param unitCell  new UnitCell<D> (with new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice type must have been set before this function is
      * called. The logical size of the FSArray<double, 6> "parameters"
      * array must match the expected number of parameters for the
      * current lattice type.
      *
      * If a space group has been declared but a basis has not yet been
      * initialized, then a symmetry-adapted basis will be constructed.
      *
      * \param parameters  array of new unit cell parameters
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      /**
      * Notify System members that unit cell parameters have been modified.
      *
      * This function should be called whenever the unit cell parameters
      * are modified. It calls functions mixture().clearUnitCellData(),
      * domain().wavelist().clearUnitCellData(), and clearCFields().
      */
      void clearUnitCellData();

      ///@}
      /// \name Field Containers
      ///@{

      /**
      * Get the monomer concentration (c) fields (const).
      */
      typename T::CFields const & c() const;

      /**
      * Get the chemical potential (w) fields (non-const).
      */
      typename T::WFields& w();

      /**
      * Get the chemical potential (w) fields (const).
      */
      typename T::WFields const & w() const;

      ///@}
      /// \name Component Object Accessors
      ///@{

      /**
      * Get the Mixture (const).
      */
      typename T::Mixture const & mixture() const;

      #if 0 // Delay implementation until needed, if ever
      /**
      * Get the MixtureModifier (non-const).
      */
      typename T::MixtureModifier& mixtureModifier();
      #endif

      /**
      * Get the %Interaction (non-const).
      */
      typename T::Interaction& interaction();

      /**
      * Get the %Interaction (const).
      */
      typename T::Interaction const & interaction() const;

      /**
      * Get the Domain (const).
      */
      typename T::Domain const & domain() const;

      #if 0 // Delay until Simulator class is finished
      /**
      * Get the Simulator (non-const).
      */
      typename T::Simulator& simulator();

      /**
      * Get the Simulator (const).
      */
      typename T::Simulator const & simulator() const;
      #endif

      /**
      * Get the FileMaster (non-const).
      *
      * Access (non-const) is used in some unit tests.
      */
      FileMaster& fileMaster();

      /**
      * Get the FileMaster (const).
      */
      FileMaster const & fileMaster() const;

      ///@}
      #if 0
      /// \name Timers
      ///@{

      /**
      * Write timer information to an output stream.
      *
      * \param out  output stream
      */
      void writeTimers(std::ostream& out) const;

      /**
      * Clear timers
      */
      void clearTimers();

      ///@}
      #endif

   protected:

      /**
      * Set the number of threads given as a command line argument.
      *
      * This function is called in the setOpts function that processes
      * command line arguments. The argument nThread may be passed to the 
      * main program as the argument of the -t option. This value gives
      * the number of threads in a threaded CPU implementation or an
      * explicit choice for the maximum number of threads per block in 
      * GPU code. 
      *
      * The do-nothing default implementation is used by CPU code that
      * has not implemented threading (the current status). 
      *
      * \param nThread  thread count
      */
      virtual void setThreadCount(int nThread)
      {};

   private:

      // Component objects

      /**
      * Chemical potential fields.
      */
      typename T::WFields w_;

      /**
      * Monomer concentration / volume fraction fields.
      */
      typename T::CFields c_;

      // Pointers to associated objects

      /**
      * Pointer to enclosing instance of System subclass.
      */
      typename T::System* systemPtr_;

      /**
      * Pointer to Mixture object (solves MDE for all species).
      */
      typename T::Mixture* mixturePtr_;

      #if 0  // Skip, for now
      /**
      * Pointer to MixtureModifier (public non-const interface for Mixture).
      */
      typename T::MixtureModifier* mixtureModifierPtr_;
      #endif

      /**
      * Pointer to Domain object (unit cell, mesh, fft, group, basis).
      */
      typename T::Domain* domainPtr_;

      /**
      * Pointer to %Interaction (excess free energy model).
      */
      typename T::Interaction* interactionPtr_;

      #if 0 // Comment out until Simulator class is finished
      /**
      * Pointer to a Simulator.
      */
      typename T::Simulator* simulatorPtr_;
      #endif

      /**
      * Filemaster (holds path prefixes for input and output files).
      */
      FileMaster* fileMasterPtr_;

      /**
      * Polymer model enumeration (thread or bead), read from file.
      */
      PolymerModel::Type polymerModel_;

      // Boolean state variables

      /**
      * Has memory been allocated for fields ?
      */
      bool isAllocated_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      // Private member functions

      /**
      * Get the Mixture by non-const reference (private interface).
      */
      typename T::Mixture & mixture_();

      /**
      * Get the Domain by non-const reference (private interface).
      */
      typename T::Domain& domain_();

      /**
      * Allocate memory for w and c fields (private).
      */
      void allocateFields();

      /**
      * Read a string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

      /**
      * Read a floating point number and echo to log file.
      *
      * Used to read numerical values in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param value  number to read and echo
      */
      void readEcho(std::istream& in, double& value) const;

   };

   // Inline member functions

   // Get the Mixture (const).
   template <int D, class T> inline 
   typename T::Mixture const & System<D,T>::mixture() const
   {  return *mixturePtr_; }

   #if 0 // Delay or skip
   // Get the MixtureModifier (non-const).
   template <int D, class T> inline 
   typename T::MixtureModifier& System<D,T>::mixtureModifier()
   {
      UTIL_ASSERT(mixtureModifierPtr_);
      return *mixtureModifierPtr_;
   }
   #endif

   // Get the %Interaction (non-const).
   template <int D, class T> inline 
   typename T::Interaction& System<D,T>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the %Interaction (const).
   template <int D, class T> inline 
   typename T::Interaction const & System<D,T>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Domain (const).
   template <int D, class T> inline 
   typename T::Domain const & System<D,T>::domain() const
   {  return *domainPtr_; }

   #if 0 // Comment out until Simulator class is finished
   // Get the Simulator (non-const).
   template <int D, class T> inline 
   typename T::Simulator& System<D,T>::simulator()
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Get the Simulator (const).
   template <int D, class T> inline 
   typename T::Simulator const & System<D,T>::simulator() const
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }
   #endif

   // Get the FileMaster (non-const).
   template <int D, class T> inline 
   FileMaster& System<D,T>::fileMaster()
   {  return *fileMasterPtr_; }

   // Get the FileMaster (const).
   template <int D, class T> inline 
   FileMaster const & System<D,T>::fileMaster() const
   {  return *fileMasterPtr_; }

   // Get the container of c fields (const).
   template <int D, class T> inline
   typename T::CFields const & System<D,T>::c() const
   {  return c_; }

   // Get the container of w fields (non-const).
   template <int D, class T> inline
   typename T::WFields& System<D,T>::w()
   {  return w_; }

   // Get the container of w fields (const).
   template <int D, class T> inline
   typename T::WFields const & System<D,T>::w() const
   {  return w_; }

   // Private inline functions:

   // Get the Mixture (non-const).
   template <int D, class T> inline 
   typename T::Mixture & System<D,T>::mixture_() 
   {  return *mixturePtr_; }

   // Get the Domain (non-const).
   template <int D, class T> inline 
   typename T::Domain & System<D,T>::domain_() 
   {  return *domainPtr_; }

} // namespace Cp
} // namespace Pscf
#endif
