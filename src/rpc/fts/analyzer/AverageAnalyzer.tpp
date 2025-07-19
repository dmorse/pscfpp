#ifndef RPC_AVERAGE_ANALYZER_TPP
#define RPC_AVERAGE_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>

#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   AverageAnalyzer<D>::AverageAnalyzer(Simulator<D>& simulator, 
                                       System<D>& system)
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&system),
      nSamplePerOutput_(1)
   {}

   /*
   * Destructor.
   */
   template <int D>
   AverageAnalyzer<D>::~AverageAnalyzer()
   {}

   /*
   * Read interval, outputFileName, and nSamplePerOutput.
   */
   template <int D>
   void AverageAnalyzer<D>::readParameters(std::istream& in)
   {
      // Read interval and outputFileName
      Analyzer<D>::readParameters(in);

      // Read nSamplePerOutput_
      nSamplePerOutput_ = 1;
      readOptional(in,"nSamplePerOutput", nSamplePerOutput_);
      if (nSamplePerOutput_ > 0) {
         std::string fileName = outputFileName(".dat");
         system().fileMaster().openOutputFile(fileName, outputFile_);
      }

      // Set the Average accumulator to compute block averages 
      // for blocks containing nSamplePerOutput_ sampled values
      accumulator_.setNSamplePerBlock(nSamplePerOutput_);
   }

   /*
   * Setup before system.
   */
   template <int D>
   void AverageAnalyzer<D>::setup()
   {  accumulator_.clear(); }

   /*
   * Compute and sample current values.
   */
   template <int D>
   void AverageAnalyzer<D>::sample(long iStep)
   {
      if (!isAtInterval(iStep)) return;

      double value = compute();
      accumulator_.sample(value);

      // Output block averages
      if (nSamplePerOutput_ > 0) {
         if (accumulator_.isBlockComplete()) {
            int beginStep = iStep - (nSamplePerOutput_ - 1)*interval();
            value = accumulator_.blockAverage();
            outputValue(beginStep, value);
         }
      }

   }

   /*
   * Write a sampled or block average value to file.
   */
   template <int D>
   void AverageAnalyzer<D>::outputValue(int step, double value) 
   {
      UTIL_CHECK(outputFile_.is_open());
      outputFile_ << Int(step);
      outputFile_ << Dbl(value);
      outputFile_ << "\n";
   }

   /*
   * Output results after a system is completed.
   */
   template <int D>
   void AverageAnalyzer<D>::output()
   {
      // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
      std::string fileName;

      #if 0
      // Write parameter (*.prm) file
      fileName = outputFileName(".prm");
      system().fileMaster().openOutputFile(fileName, outputFile_);
      ParamComposite::writeParam(outputFile_);
      outputFile_.close();
      #endif

      // Write average (*.ave) file
      fileName = outputFileName(".ave");
      system().fileMaster().openOutputFile(fileName, outputFile_);
      double ave = accumulator_.average();
      outputFile_ << "Average = ";
      outputFile_ << Dbl(ave);
      if (!simulator().hasRamp()) {
         double err = accumulator_.blockingError();
         outputFile_ << " +- " << Dbl(err, 10, 3);
      }
      outputFile_ << "\n";

      // Write error analysis to file
      if (!simulator().hasRamp()) {
         outputFile_ << "\n";
         std::string line;
         line =
         "-------------------------------------------------------------------";
         outputFile_ << line << std::endl;
         accumulator_.output(outputFile_);
         outputFile_ << std::endl;
      }
      
      outputFile_.close();

   }

}
}
#endif
