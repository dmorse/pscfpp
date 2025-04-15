
namespace Pscf{

   /**
   * \defgroup Pscf_Solver_Module Solver Templates
   *
   * Templates for "solver" classes that compute statistical mechanical 
   * properties of a gas of non-interacting molecules. 
   *
   * The templates defined in this module are designed to be used as base
   * classes for solver classes used in several different implementations
   * of polymer field theory. Different implementations may be designed 
   * for different geometries or boundary conditions or for different 
   * hardware (e.g., CPU vs. GPU), and may use different algorithms to 
   * solve the modified diffusion equation.  Source code specific to each 
   * such implementation is defined in a different implementation-level
   * subnamespace within Pscf. By convention, code that is defined in an
   * implementation level namespace may not use names or include headers 
   * for code defined in any other implementation-level namespace. This
   * convention allows the use of identical names for analogous classes 
   * in different implementation-level namespaces without causing name 
   * clashes or ambiguities.
   * 
   * The solver templates PolymerTmpl and MixtureTmpl templates are each
   * derived from a "descriptor" class that provides a description of a
   * species or mixture but that does not provide functions or data 
   * structures required to solve single-molecule statistical mechanics 
   * problem.  The descriptor classes for polymer species and mixtures are 
   * named PolymerSpecies and MixtureBase, respectively. The descriptor 
   * class for solvent species is named SolventSpecies. Source code files 
   * for these descriptor classes are in directory src/pscf/chem.
   *
   * To define an implementation of polymer field theory, one must define
   * the following set of solver classes within an implementation-level
   * subspace of Pscf:
   *
   *  - class Propagator, derived from PropagatorTmpl<Propagator>
   *  - class Block, derived from BlockTmpl<Propagator>
   *  - class Polymer, derived from PolymerTmpl<Block>
   *  - class Solvent, derived from SolventSpecies
   *  - class Mixture, derived from MixtureTmpl<Polymer, Solvent>
   *
   * \ingroup Pscf_Base_Module
   */

}
