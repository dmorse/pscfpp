#ifndef RPC_LINEAR_RESPONSE_ANALYZER_H
#define RPC_LINEAR_RESPONSE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                             // base class

#include <util/containers/DArray.h>               // member
#include <util/accumulators/Average.h>            // member
#include <prdc/cpu/RField.h>                   // member
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
   * LinearResponseAnalyzer evaluates linear reponse estimation for 
   * imcompressibility constraint. Obtain average of each step of move. 
   *
   */
   template <int D>
   class LinearResponseAnalyzer : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      LinearResponseAnalyzer(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      ~LinearResponseAnalyzer(){};

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
      * Clear accumulators.
      */
      void setup();
   
      /**
      * Add particles to BinaryStructureFactor accumulators.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      void output();
      
      /**
      * Update real error/ estimate error
      */
      void updateAccumulators();

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
      * Get Average accumulator for a specific value.
      *
      * Call only on processors that have accumulators.
      *
      * \param i integer index of value.
      */
      const Average& accumulator(int i) const;
      
      /// Number of wavevectors.
      int  nWave_;
      
      /**
      * Pointer to parent Simulator
      */
      Simulator<D>* simulatorPtr_;     
      
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_; 
      
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
      using ParamComposite::readDArray;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;

   private:

      /// Has readParam been called?
      bool isInitialized_;
      
      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      /// Fourier space dimensions
      IntVec<D> kMeshDimensions_;
      
      /// Fourier space size
      int kSize_;
      
      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;
      
      /// Field before adding random noise
      DArray< RField<D> > w1_;
      
      /// Field after adding random noise
      DArray< RField<D> > w2_;
      
      /// Actual error from incompressibility constraint in real space
      RField<D> realError_;
      
      /// Actual error from incompressibility constraint in Fourier space
      RFieldDft<D> realErrorDft_;
      
      /// Linear Response Estimation in Fourier space
      RFieldDft<D> estimateDft_;
      
      /// IntraCorrelation (homopolymer)
      RField<D> intraCorrelation_;
      
      /** 
      * Random change in pressure field with value is selected 
      * from uniform distribution [-stepSize_, stepSize_]
      */ 
      double stepSize_;
   
   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& LinearResponseAnalyzer<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& LinearResponseAnalyzer<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_LINEAR_RESPONSE_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class LinearResponseAnalyzer<1>;
   extern template class LinearResponseAnalyzer<2>;
   extern template class LinearResponseAnalyzer<3>;
   #endif

}
}
#endif
