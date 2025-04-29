#ifndef RPC_MAX_ORDER_PARAMETER_H
#define RPC_MAX_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                      // base class

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
   * MaxOrderParameter is used to detect an order-disorder transition.
   *
   * This class evalutaes maximum amplitude of the second power of the
   * Fourier mode amplitude of fluctuating fields.
   * 
   * The order parameter is defined as
   * \f[
   *     \psi(k)  = \max [ |W_{-}({\bf k})|^{2} ]
   * \f]
   * where \f$ W_{-}({\bf k})\f$ is fluctuating field component with
   * wavevector \f$ {\bf k} \f$.
   *
   * \see \ref rpc_MaxOrderParameter_page "Manual Page"
   * 
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class MaxOrderParameter : public AverageAnalyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      MaxOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      virtual ~MaxOrderParameter();
   
      /** 
      * Setup before simulation loop.
      */
      virtual void setup();
      
      using AverageAnalyzer<D>::readParameters;
      using AverageAnalyzer<D>::nSamplePerOutput;
      using AverageAnalyzer<D>::setup;
      using AverageAnalyzer<D>::sample;
      using AverageAnalyzer<D>::output;
   
   protected:

      using ParamComposite::setClassName;
      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;
      using AverageAnalyzer<D>::outputFile_;
      
      /**
      * Compute and return the max order parameter.
      */
      virtual double compute();
      
      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);
            
   private:
      
      /// Number of wavevectors in fourier mode.
      int  kSize_;
      
      /// Dimensions of fourier space 
      IntVec<D> kMeshDimensions_;
      
      /// Has readParam been called?
      bool isInitialized_;
      
      /// W_ in Fourier mode
      RFieldDft<D> wK_;

      /// Max order parameter
      double maxOrderParameter_;

   };
   
   #ifndef RPC_MAX_ORDER_PARAMETER_TPP
   // Suppress implicit instantiation
   extern template class MaxOrderParameter<1>;
   extern template class MaxOrderParameter<2>;
   extern template class MaxOrderParameter<3>;
   #endif

}
}
#endif
