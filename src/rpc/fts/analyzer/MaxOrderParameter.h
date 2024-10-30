#ifndef RPC_MAX_ORDER_PARAMETER_H
#define RPC_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                             // base class

#include <util/containers/DArray.h>               // member
#include <util/accumulators/Average.h>            // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <map>                                    // member

#include <string>
#include <iostream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * MaxOrderParameter is used to detect the ODT.
   *
   * This class evalutaes maximum amplitude of the 
   * second power of Fourier modes of fluctuating fields.
   * 
   * The order parameter is defined as
   * \f[
   *     psi(k)  = max[W_(k)W_(-k)]
   * \f]
   * where \f$W_(k)\f$ is fluctuating field at Fourier mode k 
   * 
   * \ingroup Rpc_Simulate_Analyzer_Module
   */
   template <int D>
   class MaxOrderParameter : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      ~MaxOrderParameter(){};

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);
   
      /** 
      * Setup before simulation loop.
      */
      void setup();
   
      /**
      * Compute a sampled value and update the accumulator.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      void output();
            
      /**
      * Compute max order parameter
      */
      void computeMaxOrderParameter();
   
   protected:

      /**
      * Output file stream.
      */
      std::ofstream outputFile_;
      
      /**
      * Output filename
      */
      std::string filename_;
      
      /**
      * Structure factor
      */
      double maxOrderParameter_;
   
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      /** 
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();
      
      using ParamComposite::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;

   private:

      /// Pointer to parent Simulator
      Simulator<D>* simulatorPtr_;     
      
      /// Pointer to the parent system.
      System<D>* systemPtr_; 
      
      /// Number of wavevectors in fourier mode.
      int  kSize_;
      
      /// Dimensions of fourier space 
      IntVec<D> kMeshDimensions_;
      
      /// Statistical accumulator.
      Average accumulator_;
      
      /// Whether the Average object is needed?
      bool hasAverage_;
      
      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      /// Has readParam been called?
      bool isInitialized_;
      
      /// W_ in Fourier mode
      RFieldDft<D> wK_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& MaxOrderParameter<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& MaxOrderParameter<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_MAX_ORDER_PARAMETER_TPP
   // Suppress implicit instantiation
   extern template class MaxOrderParameter<1>;
   extern template class MaxOrderParameter<2>;
   extern template class MaxOrderParameter<3>;
   #endif

}
}
#endif
