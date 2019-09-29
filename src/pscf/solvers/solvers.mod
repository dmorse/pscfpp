
namespace Pscf{

   /**
   * \defgroup Pscf_Solver_Module Solver Templates
   *
   * Templates for classes that solve modified diffusion equations.
   *
   * The templates defined in this module are designed to be used 
   * as base classes for classes that define a variety of different
   * implementations of self-consistent field theory (SCFT), in 
   * which each implementation uses a particular set of algorithms
   * to solve the modified diffusion equation (MDE) in a particular
   * type of geometry.
   *
   * To define an implementation of SCFT, one must define the 
   * following set of solver classes derived from these templates:
   *
   *  - A Propagator class, derived from PropagatorTmpl
   *  - A Block class, derived from BlockTmpl<Propagator>
   *  - A Polymer class, derived from PolymerTmpl<Block>
   *  - A Solvent class, derived from SolventTmpl<Block>
   *  - A Mixture class, derived from MixtureTmpl<Polymer, Solvent>
   *
   * \ingroup Pscf_Base_Module
   */

}
