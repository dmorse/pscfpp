#ifndef RPC_AVERAGE_ANALYZER_H
#define RPC_AVERAGE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/accumulators/Average.h>           // member

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of a single sampled real variables,
   * and optionally writes values or block averages to a data file during a 
   * simulation.  It is intended for use as a base class for any Analyzer
   * that computes and evaluates an average for a single physical variable.
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AverageAnalyzer : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Simulator object.
      * \param system parent System object.
      */
      AverageAnalyzer(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~AverageAnalyzer();

      /**
      * Read interval, outputFileName and (optionally) nSamplePerOutput.
      *
      * The optional variable nSamplePerOutput defaults to 1, which 
      * causes every sampled value to be written to file.  Setting 
      * nSamplePerOutput = 0 suppresses output of block averages to
      * file. 
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Setup before loop. 
      *
      * Opens an output file, if nSamplePerOutput > 0.
      */
      virtual void setup();

      /**
      * Compute a sampled value and update the accumulator.
      *
      * \param iStep  MD time step index
      */
      virtual void sample(long iStep);

      /**
      * Write final results to file after a simulation.
      *
      * Write an average value. If the simulation does not have a ramp,
      * it also writes an estimate of the error on the average and 
      * information estimates from hierarchichal block averages about
      * how that estimate was obtained. Information about error analysis
      * is suppressed when a ramp exists. 
      */
      virtual void output();

      /**
      * Get value of nSamplePerOutput.
      *
      * If nSamplePerOutput == 0, output of block averages is disabled.
      * For nSamplePerOutput > 0, nSamplePerOutput is the number of 
      * sampled values averaged in each block average.
      */
      int nSamplePerOutput() const;

      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;

   protected:

      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;

      /**
      * Compute value of sampled quantity.
      */
      virtual double compute() = 0;

      /**
      * Return reference to parent simulator.
      */
      Simulator<D>& simulator();

      /**
      * Return reference to parent system.
      */
      System<D>& system();

      // Output file stream
      std::ofstream outputFile_;

      /// Average object
      Average accumulator_;

   private:

      /// Pointer to the parent simulator.
      Simulator<D>* simulatorPtr_;

      /// Pointer to the parent system.
      System<D>* systemPtr_;

      /// Number of samples per block average output.
      int nSamplePerOutput_;

   };

   // Inline functions

   // Get the parent Simulator.
   template <int D>
   inline Simulator<D>& AverageAnalyzer<D>::simulator()
   {  return *simulatorPtr_; }

   // Get the parent System.
   template <int D>
   inline System<D>& AverageAnalyzer<D>::system()
   {  return *systemPtr_; }

   // Get nSamplePerOutput.
   template <int D>
   inline int AverageAnalyzer<D>::nSamplePerOutput() const
   {  return nSamplePerOutput_; }

   #ifndef RPC_AVERAGE_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class AverageAnalyzer<1>;
   extern template class AverageAnalyzer<2>;
   extern template class AverageAnalyzer<3>;
   #endif

}
}
#endif
