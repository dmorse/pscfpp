#ifndef RPG_FOURTH_ORDER_PARAMETER_H
#define RPG_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"                      // base class

#include <util/containers/DArray.h>               // member
#include <util/accumulators/Average.h>            // member
#include <prdc/cuda/RField.h>                     // member
#include <prdc/cuda/RFieldDft.h>                  // member
#include <map>                                    // member

#include <string>
#include <iostream>

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * FourthOrderParameter is used to detect an order-disorder transition.
   *
   * This class evaluates the sum of fourth power of the
   * Fourier mode amplitude of fluctuating fields.
   * 
   * The order parameter is defined as
   * \f[
   *     \Psi_{\text{fourth}} \equiv \left[ \sum W_{-}(\bf G)^4 \right] ^{\frac{1}{4}}
   * \f]
   * where \f$W_(G)\f$ is a Fourier mode of fluctuating field 
   * 
   * \ingroup Rpg_Fts_Analyzer_Module
   */
   template <int D>
   class FourthOrderParameter : public AverageAnalyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      FourthOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~FourthOrderParameter();
   
      /** 
      * Setup before simulation loop.
      */
      virtual void setup();
      
      /**
      * Compute and return the fourth order parameter.
      */
      virtual double compute();
      
      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);
            
      /**
      * Compute prefactor for each Fourier wavevector.
      * 
      * For the real-valued function W_, each Fourier 
      * coefficient G satisfies W_(G) = W_(-G). This function
      * uses Brillouin Zone (BZ) indices representation. After
      * applying fftw, if both the wavevector G and its 
      * inverse -G exist in k-space, the prefactor is
      * assigned to be 1/2 for both G and -G. Otherwise, 
      * it is assigned to be 1.
      */
      void computePrefactor();
      
      using AverageAnalyzer<D>::readParameters;
      using AverageAnalyzer<D>::nSamplePerOutput;
      using AverageAnalyzer<D>::setup;
      using AverageAnalyzer<D>::sample;
      using AverageAnalyzer<D>::output;
   
   protected:

      using AverageAnalyzer<D>::simulator;
      using AverageAnalyzer<D>::system;
      using AverageAnalyzer<D>::outputFile_;
      using ParamComposite::setClassName;
      
   private:
      
      /// Number of wavevectors in fourier mode.
      int  kSize_;
      
      /// Dimensions of fourier space 
      IntVec<D> kMeshDimensions_;
      
      /// Has readParam been called?
      bool isInitialized_;
      
      /// W_ in Fourier mode
      RFieldDft<D> wK_;
      
      /// W_ in Real space
      RField<D> wc0_;
      
      /// Prefactor for each Fourier wavevector
      RField<D> prefactor_;
      
      /// Structure factor
      double FourthOrderParameter_;

   };
   
   #ifndef RPG_FOURTH_ORDER_PARAMETER_TPP
   // Suppress implicit instantiation
   extern template class FourthOrderParameter<1>;
   extern template class FourthOrderParameter<2>;
   extern template class FourthOrderParameter<3>;
   #endif

}
}
#endif
