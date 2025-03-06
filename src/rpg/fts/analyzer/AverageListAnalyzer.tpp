#ifndef RPG_AVERAGE_LIST_ANALYZER_TPP
#define RPG_AVERAGE_LIST_ANALYZER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"

#include <rpg/System.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

namespace Pscf {
namespace Rpg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   AverageListAnalyzer<D>::AverageListAnalyzer(System<D>& system)
    : Analyzer<D>(),
      systemPtr_(&system),
      nSamplePerOutput_(1),
      nValue_(0),
      hasAccumulators_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   AverageListAnalyzer<D>::~AverageListAnalyzer()
   {}

   /*
   * Read interval and outputFileName.
   */
   template <int D>
   void AverageListAnalyzer<D>::readParameters(std::istream& in)
   {
      Analyzer<D>::readParameters(in);
      nSamplePerOutput_ = 1;
      readOptional(in,"nSamplePerOutput", nSamplePerOutput_);
      if (nSamplePerOutput() > 0) {
         std::string fileName = outputFileName(".dat");
         system().fileMaster().openOutputFile(fileName, outputFile_);
      }
      // Note: ReadParameters method of derived classes should call this,
      // determine nValue and then call initializeAccumulators(nValue).
   }

   /*
   * Clear accumulators (do nothing on slave processors).
   */
   template <int D>
   void AverageListAnalyzer<D>::clear()
   {
      UTIL_CHECK(hasAccumulators_);
      clearAccumulators();
   }

   /*
   * Setup before system.
   */
   template <int D>
   void AverageListAnalyzer<D>::setup()
   {
      UTIL_CHECK(hasAccumulators_);
      clear();
   }

   /*
   * Compute and sample current values.
   */
   template <int D>
   void AverageListAnalyzer<D>::sample(long iStep)
   {
      UTIL_CHECK(hasAccumulators_);
      if (!isAtInterval(iStep)) return;
      compute();
      updateAccumulators(iStep);
   }

   /*
   * Output results after a system is completed.
   */
   template <int D>
   void AverageListAnalyzer<D>::output()
   {
      UTIL_CHECK(hasAccumulators_);

      // Close data file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      #if 0
      // Write parameter (*.prm) file
      system().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_);
      outputFile_.close();
      #endif

      // Write average (*.ave) and error analysis (*.aer) files
      outputAccumulators();

   }

   /**
   * Set nValue and allocate arrays with dimensions nValue.
   */
   template <int D>
   void AverageListAnalyzer<D>::initializeAccumulators(int nValue)
   {
      UTIL_CHECK(nValue > 0);
      UTIL_CHECK(nValue_ == 0);
      UTIL_CHECK(nSamplePerOutput_ >= 0);

      // Allocate arrays
      accumulators_.allocate(nValue);
      names_.allocate(nValue);
      values_.allocate(nValue);
      nValue_ = nValue;
      hasAccumulators_ = true;

      // Set the accumulators to compute block averages with
      // nSamplePerOutput_ sampled values per block
      for (int i = 0; i < nValue_; ++i) {
         accumulators_[i].setNSamplePerBlock(nSamplePerOutput_);
      }

      clearAccumulators();
   }

   /*
   * Clear accumulators.
   */
   template <int D>
   void AverageListAnalyzer<D>::clearAccumulators()
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(nValue_ > 0);
      for (int i = 0; i < nValue_; ++i) {
         accumulators_[i].clear();
      }
   }

   template <int D>
   void AverageListAnalyzer<D>::setName(int i, std::string name)
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      names_[i] = name;
   }

   /*
   * Update accumulators for all current values.
   */
   template <int D>
   void AverageListAnalyzer<D>::updateAccumulators(long iStep)
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(accumulators_.capacity() == nValue_);

      // Update accumulators.
      for (int i = 0; i < nValue(); ++i) {
         double data = value(i);
         accumulators_[i].sample(data);
      }

      // Output block averages
      if (nSamplePerOutput_ > 0) {
         if (accumulators_[0].isBlockComplete()) {
            UTIL_CHECK(outputFile_.is_open());
            int beginStep = iStep - (nSamplePerOutput_ - 1)*interval();
            outputFile_ << Int(beginStep);
            for (int i = 0; i < nValue(); ++i) {
               UTIL_CHECK(accumulators_[i].isBlockComplete());
               double block = accumulators_[i].blockAverage();
               outputFile_ << Dbl(block);
            }
            outputFile_ << "\n";
         }
      }

   }

   /*
   * Output results to file after simulation is completed.
   */
   template <int D>
   void AverageListAnalyzer<D>::outputAccumulators()
   {
      UTIL_CHECK(hasAccumulators_);

      // Close data (*.dat) file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      // Compute maximum length of name field
      int nameWidth = 0;
      int length;
      for (int i = 0; i < nValue_; ++i) {
         length = names_[i].length();
         if (length > nameWidth) {
            nameWidth = length;
         }
      }
      nameWidth += 2;

      // Write average (*.ave) file
      std::string fileName = outputFileName(".ave");
      system().fileMaster().openOutputFile(fileName, outputFile_);
      double ave, err;
      for (int i = 0; i < nValue_; ++i) {
         ave = accumulators_[i].average();
         err = accumulators_[i].blockingError();
         outputFile_ << " " << std::left << std::setw(nameWidth)
                      << names_[i] << "   ";
         outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      outputFile_.close();

      // Write error analysis (*.aer) file
      fileName = outputFileName(".aer");
      system().fileMaster().openOutputFile(fileName, outputFile_);
      std::string line;
      line =
      "---------------------------------------------------------------------";
      for (int i = 0; i < nValue_; ++i) {
         outputFile_ << line << std::endl;
         outputFile_ << names_[i] << " :" << std::endl;
         accumulators_[i].output(outputFile_);
         outputFile_ << std::endl;
      }
      outputFile_.close();

      #if 0
      // Write data format file (*.dfm) file
      fileName = outputFileName();
      fileName += ".dfm";
      system().fileMaster().openOutputFile(fileName, outputFile_);
      outputFile_ << "Value = " << nValue() << std::endl;
      outputFile_ << "iStep  ";
      for (int i = 0; i < nValue_; ++i) {
         outputFile_ << names_[i] << "  ";
      }
      outputFile_ << std::endl;
      outputFile_.close();
      #endif

   }

}
}
#endif
