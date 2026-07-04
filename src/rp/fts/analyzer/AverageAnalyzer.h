#ifndef RP_AVERAGE_ANALYZER_H
#define RP_AVERAGE_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/Analyzer.h>     // base class
#include <util/accumulators/Average.h>    // member
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Analyze averages and block averages of a single real variable.
   *
   * This class evaluates the average of a single sampled real variables,
   * and optionally writes values or block averages to a data file during a
   * simulation.  It is intended for use as a base class for any Analyzer
   * that computes and evaluates an average for a single physical variable.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>.
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class AverageAnalyzer : public Analyzer<D,T>
   {

   public:

      // Protected constructor and destructor (see below).

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

   protected:

      // Protected member functions

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      AverageAnalyzer(Simulator<D,T>& simulator, 
                      System<D,T>& system);

      /**
      * Destructor.
      */
      ~AverageAnalyzer() = default;

      /**
      * Compute value of sampled quantity.
      */
      virtual double compute() = 0;

      /**
      * Output a sampled or block average value.
      *
      * \param step  value for step counter
      * \param value  value of physical observable
      */
      virtual void outputValue(int step, double value);

      // Protected member variables

      /// Output file stream.
      std::ofstream outputFile_;

      /// Average object.
      Average accumulator_;

   private:

      /// Number of samples per block average output.
      int nSamplePerOutput_;

      using AnalyzerT = Analyzer<D,T>;

   };

   // Inline functions

   // Get nSamplePerOutput.
   template <int D, class T> inline 
   int AverageAnalyzer<D,T>::nSamplePerOutput() const
   {  return nSamplePerOutput_; }

}
}
#endif
