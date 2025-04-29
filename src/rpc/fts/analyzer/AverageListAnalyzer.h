#ifndef RPC_AVERAGE_LIST_ANALYZER_H
#define RPC_AVERAGE_LIST_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <util/accumulators/Average.h>           // member

namespace Pscf {
namespace Rpc
{

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of several sampled real variables, and
   * optionally writes block averages to a data file during a simulation.
   * It is intended for use as a base class for Analyzers that evaluate
   * averages and (optionally) block averages for several physical variables.
   *
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class AverageListAnalyzer : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent simulator object.
      * \param system  parent system object.
      */
      AverageListAnalyzer(Simulator<D>& simulator, System<D>& system);

      /**
      * Destructor.
      */
      virtual ~AverageListAnalyzer();

      /**
      * Read interval, outputFileName and (optionally) nSamplePerOutput.
      *
      * The optional variable nSamplePerOutput defaults to 0, which 
      * disables computation and output of block averages. Setting 
      * nSamplePerOutput = 1 outputs every sampled value.
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
      */
      virtual void output();

      /**
      * Get value of nSamplePerOutput.
      *
      * If nSamplePerOutput == 0, output of block averages is disabled.
      * For nSamplePerOutput > 0, the value is the number of sampled values
      * averaged in each block.
      */
      int nSamplePerOutput() const;

      /**
      * Get number of variables.
      */
      int nValue() const;

      /**
      * Get name associated with value.
      *
      * Call only on processors that have accumulators.
      *
      * \param i integer index of name/value pair.
      */
      const std::string& name(int i) const;

      /**
      * Get Average accumulator for a specific value.
      *
      * \param i integer index of value.
      */
      const Average& accumulator(int i) const;

      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;

   protected:

      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;

      /**
      * Set name of variable.
      *
      * Call only on master.
      *
      * \param i integer index of variable
      * \param name name of variable number i
      */
      void setName(int i, std::string name);

      /**
      * Instantiate Average accumulators and set nValue.
      *
      * \pre hasAccumulator == false
      * \pre nSamplePerOutput >= 0
      */
      void initializeAccumulators(int nValue);

      /**
      * Clear internal state of all accumulators.
      *
      * \pre hasAccumulator == true
      */
      void clearAccumulators();

      /**
      * Compute value of sampled quantities.
      *
      * Call on all processors.
      */
      virtual void compute() = 0;

      /**
      * Set current value, used by compute function.
      *
      * \param i integer index of variable
      * \param value current value of variable
      */
      void setValue(int i, double value);

      /**
      * Get current value of a specific variable.
      *
      * Call only on master.
      *
      * \param i integer index of variable.
      */
      double value(int i) const;

      /**
      * Add current value to accumulator, output block average if needed.
      *
      * \param iStep simulation step counter
      */
      void updateAccumulators(long iStep);

      /**
      * Write results of statistical analysis to files.
      */
      void outputAccumulators();

      /** 
      * Return reference to parent Simulator.
      */
      Simulator<D>& simulator();

      /**
      * Return reference to parent system.
      */
      System<D>& system();

      /// Output file stream
      std::ofstream outputFile_;

   private:

      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;

      /// Array of current values (only allocated on master processor)
      DArray<double> values_;

      /// Array of value names (only allocated on master processor)
      DArray<std::string> names_;

      /// Pointer to parent Simulator
      Simulator<D>* simulatorPtr_;     

      /// Pointer to the parent system.
      System<D>* systemPtr_;
 
      /// Number of samples per block average output.
      int nSamplePerOutput_;

      /// Number of values.
      int nValue_;

      /// Does this processor have accumulators ?
      bool hasAccumulators_;

   };

   // Inline functions

   /*
   * Get nSamplePerOutput.
   */
   template <int D>
   inline int AverageListAnalyzer<D>::nSamplePerOutput() const
   {  return nSamplePerOutput_; }

   /*
   * Get nValue (number of variables).
   */
   template <int D>
   inline int AverageListAnalyzer<D>::nValue() const
   {
      UTIL_CHECK(hasAccumulators_);
      return nValue_;
   }

   /*
   * Get current value of a variable, set by compute function.
   */
   template <int D>
   inline double AverageListAnalyzer<D>::value(int i) const
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      return values_[i];
   }

   /*
   * Get name of specific variable.
   */
   template <int D>
   inline const std::string& AverageListAnalyzer<D>::name(int i) const
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      return names_[i];
   }

   /*
   * Get accumulator associated with a variable.
   */
   template <int D>
   inline const Average& AverageListAnalyzer<D>::accumulator(int i) const
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      return accumulators_[i];
   }

   /*
   * Set current value of a variable.
   */
   template <int D>
   inline void AverageListAnalyzer<D>::setValue(int i, double value)
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      values_[i] = value;
   }

   // Get the parent system.
   template <int D>
   inline System<D>& AverageListAnalyzer<D>::system()
   {  return *systemPtr_; }
   
   // Get parent Simulator object.
   template <int D>
   inline Simulator<D>& AverageListAnalyzer<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_AVERAGE_LIST_ANALYZER_TPP
   // Suppress implicit instantiation
   extern template class AverageListAnalyzer<1>;
   extern template class AverageListAnalyzer<2>;
   extern template class AverageListAnalyzer<3>;
   #endif

}
}
#endif
