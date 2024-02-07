
namespace Pscf{

   /**
   * \defgroup Pscf_Solver_Module Solver Templates
   *
   * Templates for classes that solve modified diffusion equations.
   *
   * The templates defined in this module are designed to be used 
   * as base classes for classes that define a variety of different
   * implementations of polymer field theory. Different implementations
   * may be designed for different geometries or boundary conditions
   * or for different hardware (e.g., CPU vs. GPU), and may use 
   * different algorithms to solve the modified diffusion equation.  
   * Source code specific to each such implementation is defined in 
   * a different enclosed namespace of Pscf, thus allowing the use 
   * of identical names for analogous classes in different 
   * implementations without causing name clashes or ambiguities.
   *
   * To define an implementation of polymer field theory, one must 
   * define the following set of solver classes derived from these 
   * solver templates:
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
