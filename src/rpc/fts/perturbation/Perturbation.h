#ifndef RPC_PERTURBATION_H
#define RPC_PERTURBATION_H

#include <util/param/ParamComposite.h>      // base class

namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Prdc{
namespace Cpu{
   template <int D> class RField;
}
}
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   template <int D> class System;
   template <int D> class Simulator;

   /**
   * Base class for additive perturbations of standard FTS Hamiltonian.
   *
   * \ingroup Rpc_Simulate_Perturbation_Module
   */
   template <int D>
   class Perturbation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Perturbation(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~Perturbation();

      /**
      * Read parameters from archive.
      *
      * Empty default implementation.
      *
      * \param in input parameter stream
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
      * \param unperturbedHamiltonian Hamiltonian without perturbation
      */
      virtual double hamiltonian(double unperturbedHamiltonian);

      /**
      * Modify the generalized forces to include perturbation.
      *
      * Empty default implementation.
      */
      virtual void incrementDc(DArray< RField<D> >& dc);

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
      * Get parent Simulator<D> by const reference.
      */
      Simulator<D> const & simulator() const;
      
      /** 
      * Get parent System<D> by const reference.
      */      
      System<D> const & system() const;

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
      * \param lambda  new value for lambda perturbation parameter
      */
      void setLambda(double lambda);

   protected:

      /**
      * Get parent Simulator<D> by non-const reference.
      */
      Simulator<D>& simulator();
      
      /**
      * Get parent system<D> by non-const reference.
      */
      System<D>& system();
      
      /**
      * Strength of the perturbation
      */
      double lambda_;

   private:

      /// Pointer to parent Simulator.
      Simulator<D>* simulatorPtr_;
      
      /// Pointer to parent System.
      System<D>* systemPtr_;

   };

   // Inline methods

   // Return parent simulator by const reference.
   template <int D>
   inline Simulator<D> const & Perturbation<D>::simulator() const
   {
      assert(simulatorPtr_);  
      return *simulatorPtr_; 
   }

   // Return parent simulator by non-const reference.
   template <int D>
   inline Simulator<D> & Perturbation<D>::simulator() 
   {  
      assert(simulatorPtr_);  
      return *simulatorPtr_; 
   }
   
   // Return parent simulator by const reference.
   template <int D>
   inline System<D> const & Perturbation<D>::system() const
   {  
      assert(systemPtr_);  
      return *systemPtr_; 
   }
   
   // Return parent simulator by non-const reference.
   template <int D>
   inline System<D> & Perturbation<D>::system() 
   {  
      assert(systemPtr_);  
      return *systemPtr_; 
   }

   #ifndef RPC_PERTURBATION_TPP
   // Suppress implicit instantiation
   extern template class Perturbation<1>;
   extern template class Perturbation<2>;
   extern template class Perturbation<3>;
   #endif

}
}
#endif
