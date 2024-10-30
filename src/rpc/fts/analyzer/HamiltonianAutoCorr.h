#ifndef RPC_HAMILTONIAN_AUTO_CORR_H
#define RPC_HAMILTONIAN_AUTO_CORR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/
#include "Analyzer.h"                             // base class
#include <util/accumulators/AutoCorrelation.h>  // member template 
#include <util/accumulators/Average.h>            // member template

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Compute hamiltonian autocorrelation.
   *
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class HamiltonianAutoCorr : public Analyzer<D>
   {

   public:
   
      /**
      * Constructor.
      */
      HamiltonianAutoCorr(Simulator<D>& simulator, System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~HamiltonianAutoCorr()
      {} 
      
      /** 
      * Setup before beginning of loop.
      */
      void setup();
   
      /**
      * Sample hamiltonian and add to accumulator.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      void output();
      
      /**
      * Read interval and output file name.
      *
      * \param in  input parameter file
      */
      void readParameters(std::istream& in);
      
      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;     
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_; 
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;
      
   protected:      
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /** 
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();
 
   private: 
      /// Output file stream
      std::ofstream outputFile_;

      /// Output filename
      std::string filename_;

      /// Statistical accumulator.
      AutoCorrelation<double, double>  accumulator_;
      Average avgAccumulator_;
 
      /// Number of samples per block average output (number of values stored in buffer)
      int  bufferCapacity_;

      /// Maximum id for AutoCorrStage
      int  maxStageId_;

      /// BlockFactor for AutoCorrStage algorithm
      int  blockFactor_;
      
      /// Has readParam been called?
      bool isInitialized_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& HamiltonianAutoCorr<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& HamiltonianAutoCorr<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_HAMILTONIAN_AUTO_CORR_TPP
   // Suppress implicit instantiation
   extern template class HamiltonianAutoCorr<1>;
   extern template class HamiltonianAutoCorr<2>;
   extern template class HamiltonianAutoCorr<3>;
   #endif

}
}
#endif 
