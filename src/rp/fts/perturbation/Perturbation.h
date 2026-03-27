#ifndef RP_PERTURBATION_H
#define RP_PERTURBATION_H

#include <util/param/ParamComposite.h>      // base class

// Forward declaration
namespace Util {  template <typename T> class DArray; }

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Base class for additive perturbations of standard FTS Hamiltonian.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named Perturbation, that
   * are defined in Rpc and Rpg namespaces for use in the pscf_rpc and
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref psfts_perturb_page "Manual Page"
   * \ingroup Rp_Fts_Perturbation_Module
   */
   template <int D, class T>
   class Perturbation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      Perturbation(typename T::Simulator& simulator);

      /**
      * Destructor.
      */
      virtual ~Perturbation();

      /**
      * Read parameters from archive.
      *
      * Empty default implementation.
      *
      * \param in  input parameter file stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Complete any required initialization.
      *
      * This method must be called just before the beginning of
      * the main simulation loop, after an initial configuration
      * is known. It may be used to complete any initialization
      * that cannot be completed in the readParameters function.
      *
      * The default implementation is an empty function.
      */
      virtual void setup();

      /**
      * Compute and return the perturbation to the Hamiltonian.
      *
      * Default implementation returns 0.
      *
      * \param unperturbedHamiltonian  Hamiltonian without perturbation
      */
      virtual double hamiltonian(double unperturbedHamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      *
      * Empty default implementation.
      */
      virtual void incrementDc(DArray< typename T::RField >& dc);

      /**
      * Save any required internal state variables.
      *
      * This function should save any state variables that would need to
      * be restored after a rejected Monte Carlo move or failure of the
      * compressor to converge after an attempted Brownian dynamics move.
      */
      virtual void saveState();

      /**
      * Restore any required internal state variables.
      *
      * This function is called after rejection of an MC move or failure
      * of an attempted BD step, and should restore the variables saved
      * by the saveState function.
      */
      virtual void restoreState();

      /**
      * Compute and return derivative of H w/ respect to parameter lambda.
      *
      * Default implementation returns 0.
      */
      virtual double df();

      /**
      * Get parent typename T::Simulator by const reference.
      */
      typename T::Simulator const & simulator() const;

      /**
      * Get parent typename T::System by const reference.
      */
      typename T::System const & system() const;

      /**
      * Get the perturbation parameter.
      *
      * The perturbation parameter lambda is initialized to 1.0 in
      * the Perturbation constructor.
      */
      double lambda() const
      {  return lambda_; }

      /**
      * Set the perturbation parameter value.
      *
      * \param lambda  new value for the perturbation parameter
      */
      void setLambda(double lambda);

   protected:

      /**
      * Get parent Simulator by non-const reference.
      */
      typename T::Simulator& simulator();

      /**
      * Get parent System by non-const reference.
      */
      typename T::System& system();

      /**
      * Strength of the perturbation
      */
      double lambda_;

   private:

      /// Pointer to parent Simulator.
      typename T::Simulator* simulatorPtr_;

      /// Pointer to parent System.
      typename T::System* systemPtr_;

   };

   // Inline methods

   // Return parent simulator by const reference.
   template <int D, class T>
   inline typename T::Simulator const & Perturbation<D,T>::simulator() const
   {
      assert(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Return parent simulator by non-const reference.
   template <int D, class T>
   inline typename T::Simulator & Perturbation<D,T>::simulator()
   {
      assert(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Return parent simulator by const reference.
   template <int D, class T>
   inline typename T::System const & Perturbation<D,T>::system() const
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

   // Return parent simulator by non-const reference.
   template <int D, class T>
   inline typename T::System & Perturbation<D,T>::system()
   {
      assert(systemPtr_);
      return *systemPtr_;
   }

}
}
#endif
