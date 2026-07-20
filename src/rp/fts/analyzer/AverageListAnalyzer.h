#ifndef RP_AVERAGE_LIST_ANALYZER_H
#define RP_AVERAGE_LIST_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/Analyzer.h>    // base class
#include <util/accumulators/Average.h>   // member
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class template evaluates averages of several sampled real 
   * variables, and optionally writes block averages to a data file during 
   * a simulation. It is intended for use as a base class for Analyzers
   * that evaluate averages and (optionally) block averages for several
   * physical variables.
   *
   * Template parameters:
   *
   *    - D : dimension of space
   *    - T : collection of type aliases, e.g., CppTp<D>
   *
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class AverageListAnalyzer : public Analyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      AverageListAnalyzer(Simulator<D,T>& simulator, 
                          System<D,T>& system);

      /**
      * Read interval, outputFileName and (optionally) nSamplePerOutput.
      *
      * The optional variable nSamplePerOutput defaults to 0, which
      * disables computation and output of block averages. Setting
      * nSamplePerOutput = 1 outputs every sampled value.
      *
      * \param in  input parameter stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Setup before loop.
      *
      * Opens an output file, if nSamplePerOutput > 0.
      */
      void setup() override;

      /**
      * Compute sampled values and update the accumulators.
      *
      * \param iStep  simulation step index
      */
      void sample(long iStep) override;

      /**
      * Write final results to file after a simulation.
      */
      void output() override;

      /**
      * Get value of nSamplePerOutput.
      *
      * If nSamplePerOutput == 0, output of block averages is disabled.
      * For nSamplePerOutput > 0, the value is the number of sampled
      * values averaged in each block.
      */
      int nSamplePerOutput() const;

      /**
      * Get the number of variables.
      */
      int nValue() const;

      /**
      * Get name associated with a variable.
      *
      * \param i  integer index of name/value pair
      */
      const std::string& name(int i) const;

      /**
      * Get Average accumulator for a specific variable.
      *
      * \param i  integer index of value
      */
      const Average& accumulator(int i) const;

   protected:

      /// Output file stream.
      std::ofstream outputFile_;

      /**
      * Initialize Average accumulators and set nValue.
      *
      * \pre hasAccumulator == false
      * \pre nSamplePerOutput >= 0
      *
      * \param nValue  number of values
      */
      void initializeAccumulators(int nValue);

      /**
      * Clear internal state of all accumulators.
      *
      * \pre hasAccumulator == true
      */
      void clearAccumulators();

      /**
      * Set name of variable.
      *
      * \param i  integer index of variable
      * \param name  name of variable number i
      */
      void setName(int i, std::string name);

      /**
      * Set current value, used by compute function.
      *
      * \param i  integer index of variable
      * \param value  current value of variable
      */
      void setValue(int i, double value);

      /**
      * Compute values of sampled quantities.
      */
      virtual void compute() = 0;

      /**
      * Get current value of a specific variable.
      *
      * \param i  integer index of variable
      */
      double value(int i) const;

      /**
      * Add current value to accumulator, output block average if needed.
      *
      * \param iStep  simulation step counter
      */
      void updateAccumulators(long iStep);

      /**
      * Write results of statistical analysis to files.
      */
      void outputAccumulators();

      using AnalyzerT = Analyzer<D,T>;
      using AnalyzerT::simulator;
      using AnalyzerT::system;

   private:

      /// Array of Average accumulator objects.
      DArray<Average> accumulators_;

      /// Array of current values.
      DArray<double> values_;

      /// Array of variable names.
      DArray<std::string> names_;

      /// Number of samples per block average output.
      int nSamplePerOutput_;

      /// Number of variables.
      int nValue_;

      /// Does this processor have accumulators ?
      bool hasAccumulators_;

   };

   // Inline functions

   /*
   * Get nSamplePerOutput.
   */
   template <int D, class T> inline
   int AverageListAnalyzer<D,T>::nSamplePerOutput() const
   {  return nSamplePerOutput_; }

   /*
   * Get nValue (number of variables).
   */
   template <int D, class T> inline
   int AverageListAnalyzer<D,T>::nValue() const
   {
      UTIL_CHECK(hasAccumulators_);
      return nValue_;
   }

   /*
   * Get current value of a variable, set by compute function.
   */
   template <int D, class T> inline
   double AverageListAnalyzer<D,T>::value(int i) const
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      return values_[i];
   }

   /*
   * Get name of specific variable.
   */
   template <int D, class T> inline
   const std::string& AverageListAnalyzer<D,T>::name(int i) const
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      return names_[i];
   }

   /*
   * Get accumulator associated with a variable.
   */
   template <int D, class T> inline
   const Average& AverageListAnalyzer<D,T>::accumulator(int i) const
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      return accumulators_[i];
   }

   /*
   * Set current value of a variable.
   */
   template <int D, class T> inline
   void AverageListAnalyzer<D,T>::setValue(int i, double value)
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      values_[i] = value;
   }

}
}
#endif
