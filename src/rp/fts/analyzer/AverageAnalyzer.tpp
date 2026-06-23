#ifndef RP_AVERAGE_ANALYZER_TPP
#define RP_AVERAGE_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"

#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   AverageAnalyzer<D,T>::AverageAnalyzer(Simulator<D,T>& simulator,
                                         System<D,T>& system)
    : AnalyzerT(simulator, system),
      nSamplePerOutput_(1)
   {  AnalyzerT::setFileMaster(system.fileMaster()); }

   /*
   * Read interval, outputFileName, and nSamplePerOutput.
   */
   template <int D, class T>
   void AverageAnalyzer<D,T>::readParameters(std::istream& in)
   {
      // Read interval and outputFileName
      AnalyzerT::readParameters(in);

      // Read nSamplePerOutput_
      nSamplePerOutput_ = 1;
      ParamComposite::readOptional(in,"nSamplePerOutput", 
                                   nSamplePerOutput_);
      if (nSamplePerOutput_ > 0) {
         std::string fileName = AnalyzerT::outputFileName(".dat");
         AnalyzerT::system().fileMaster().openOutputFile(fileName, outputFile_);
      }

      // Set the Average accumulator to compute block averages
      // for blocks containing nSamplePerOutput_ sampled values
      accumulator_.setNSamplePerBlock(nSamplePerOutput_);
   }

   /*
   * Setup before system.
   */
   template <int D, class T>
   void AverageAnalyzer<D,T>::setup()
   {  accumulator_.clear(); }

   /*
   * Compute and sample current values.
   */
   template <int D, class T>
   void AverageAnalyzer<D,T>::sample(long iStep)
   {
      if (!AnalyzerT::isAtInterval(iStep)) return;

      double value = compute();
      accumulator_.sample(value);

      // Output value or block average
      if (nSamplePerOutput_ > 0) {
         if (nSamplePerOutput_ == 1) {
            outputValue(iStep, value);
         } else 
         if (accumulator_.isBlockComplete()) {
            int interval = AnalyzerT::interval();
            int beginStep = iStep - (nSamplePerOutput_ - 1)*interval;
            value = accumulator_.blockAverage();
            outputValue(beginStep, value);
         }
      }

   }

   /*
   * Write a sampled or block average value to file.
   */
   template <int D, class T>
   void AverageAnalyzer<D,T>::outputValue(int step, double value)
   {
      UTIL_CHECK(outputFile_.is_open());
      outputFile_ << Int(step);
      outputFile_ << Dbl(value);
      outputFile_ << "\n";
   }

   /*
   * Output results after a system is completed.
   */
   template <int D, class T>
   void AverageAnalyzer<D,T>::output()
   {
      // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }
      std::string fileName;

      #if 0
      // Write parameter (*.prm) file
      fileName = AnalyzerT::outputFileName(".prm");
      AnalyzerT::system().fileMaster().openOutputFile(fileName, 
                                                      outputFile_);
      ParamComposite::writeParam(outputFile_);
      outputFile_.close();
      #endif

      // Write average (*.ave) file
      fileName = AnalyzerT::outputFileName(".ave");
      AnalyzerT::system().fileMaster().openOutputFile(fileName, 
                                                      outputFile_);
      double ave = accumulator_.average();
      outputFile_ << "Average = ";
      outputFile_ << Dbl(ave);
      if (!AnalyzerT::simulator().hasRamp()) {
         double err = accumulator_.blockingError();
         outputFile_ << " +- " << Dbl(err, 10, 3);
      }
      outputFile_ << "\n";

      // Write error analysis to file
      if (!AnalyzerT::simulator().hasRamp()) {
         outputFile_ << "\n";
         std::string line;
         line =
         "----------------------------------------------------------------";
         outputFile_ << line << std::endl;
         accumulator_.output(outputFile_);
         outputFile_ << std::endl;
      }

      outputFile_.close();
   }

} // namespace Rp
} // namespace Rpc
#endif
