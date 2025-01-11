#ifndef RPG_FOURTH_ORDER_PARAMETER_H
#define RPG_FOURTH_ORDER_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                             // base class

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
   * FourthOrderParameter is used to detect the ODT.
   *
   * This class evalutaes the sum of fourth power of 
   * fluctuating field in fourier mode
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
   class FourthOrderParameter : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      FourthOrderParameter(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      ~FourthOrderParameter(){};

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
      * Compute fourth order parameter
      */
      void computeFourthOrderParameter();
      
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
      double FourthOrderParameter_;
   
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
      
      /// W_ in Real space
      RField<D> wc0_;
      
      /// Prefactor for each Fourier wavevector
      RField<D> prefactor_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& FourthOrderParameter<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& FourthOrderParameter<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPG_FOURTH_ORDER_PARAMETER_TPP
   // Suppress implicit instantiation
   extern template class FourthOrderParameter<1>;
   extern template class FourthOrderParameter<2>;
   extern template class FourthOrderParameter<3>;
   #endif

}
}
#endif
