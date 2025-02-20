#ifndef RPC_THERMO_DERIVATIVE_ANALYZER_TPP
#define RPC_THERMO_DERIVATIVE_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ThermoDerivativeAnalyzer.h"

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/ramp/Ramp.h>

#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <string>
#include <algorithm>

namespace Pscf {
namespace Rpc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ThermoDerivativeAnalyzer<D>::ThermoDerivativeAnalyzer(
                                Simulator<D>& simulator, System<D>& system)
    : Analyzer<D>(),
      outputFileName_(""),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      hasAverage_(true),
      hasOutputFile_(false)
   { setClassName("ThermoDerivativeAnalyzer"); }

   /*
   * Destructor.
   */
   template <int D>
   ThermoDerivativeAnalyzer<D>::~ThermoDerivativeAnalyzer()
   {}

   /*
   * Read interval and outputFileName.
   */
   template <int D>
   void ThermoDerivativeAnalyzer<D>::readParameters(std::istream& in)
   {
      readInterval(in);
      read(in, "outputFileName", outputFileName_);
      readOptional(in, "hasAverage", hasAverage_);

      if (!outputFileName_.empty()) {
         hasOutputFile_ = true;
      }

      if (hasOutputFile_){
         system().fileMaster().openOutputFile(outputFileName_, outputFile_);
         outputFile_ << "    Variable       " << "Derivative" << "\n";

      }
   }

   /*
   * Setup before simulation loop.
   */
   template <int D>
   void ThermoDerivativeAnalyzer<D>::setup()
   {
      if (hasAverage_){
         accumulator_.clear();
      }

   }

   template <int D>
   void ThermoDerivativeAnalyzer<D>::sample(long iStep)
   {
      if (!isAtInterval(iStep)) return;

      double derivative;
      derivative = computeDerivative();
      if (hasAverage_){
         accumulator_.sample(derivative);
      }

      if (hasOutputFile_){
          UTIL_CHECK(outputFile_.is_open());
          outputFile_ << Dbl(variable());
          outputFile_ << Dbl(derivative);
          outputFile_<< "\n";
      }
   }

   /*
   * Output results after a simulation is completed.
   */
   template <int D>
   void ThermoDerivativeAnalyzer<D>::output()
   {
      if (hasAverage_){
         Log::file() << std::endl;
         Log::file() << parameterType() << " : "<< std::endl;
         Log::file() << "Average of derivative : "
                     << Dbl(accumulator_.average())
                     << " +- " << Dbl(accumulator_.blockingError(), 9, 2)
                     << "\n";
      }

      // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      // Write average (*.ave) file and error analysis (*.aer) file
      if (hasAverage_){

         // Output average results to average file
         std::string fileName;
         std::string type;
         fileName = outputFileName_+  ".ave";
         system().fileMaster().openOutputFile(fileName, outputFile_);
         double ave, err;
         ave = accumulator_.average();
         err = accumulator_.blockingError();
         outputFile_ << parameterType() << "   ";
         outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         outputFile_.close();

         // Write error analysis (*.aer) file, if there is no Ramp
         if (!simulator().hasRamp()) {
            fileName = outputFileName_+  ".aer";
            system().fileMaster().openOutputFile(fileName, outputFile_);
            outputFile_ << parameterType() << ":"<< std::endl;
            accumulator_.output(outputFile_);
            outputFile_ << std::endl;
            outputFile_.close();
         }

      }

   }

}
}
#endif
