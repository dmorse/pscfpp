#ifndef PSPG_AVERAGE_LIST_ANALYZER_H
#define PSPG_AVERAGE_LIST_ANALYZER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"
#include <rpg/System.h>
#include <util/accumulators/Average.h>           // member

namespace Pscf {
namespace Rpg
{

   template <int D> class System;

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of several sampled real variables, and
   * optionally writes block averages to a data file during a simulation. 
   * It is intended for use as a base class for Analyzers that evaluate 
   * averages and (optionally) block averages for specific physical variables.
   *
   * \ingroup Rpg_Simulate_Analyzer_Module
   */
   template <int D>
   class AverageListAnalyzer : public Analyzer<D>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent SystemType object. 
      */
      AverageListAnalyzer(System<D>& system);
   
      /**
      * Destructor.
      */
      virtual ~AverageListAnalyzer(); 

      /**
      * Read interval, outputFileName and (optionally) nSamplePerBlock.
      *
      * The optional variable nSamplePerBlock defaults to 0, which disables
      * computation and output of block averages. Setting nSamplePerBlock = 1
      * outputs every sampled value. 
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      #if 0
      /**
      * Load internal state from an input archive.
      *
      * \param ar  input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an output archive.
      *
      * \param ar  output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
      #endif
      
      /**
      * Clear accumulator on master, do nothing on other processors.
      */
      virtual void clear();
  
      /**
      * Setup before loop. Opens an output file, if any.
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
      * Get value of nSamplePerBlock.
      * 
      * If nSamplePerBlock == 0, output of block averages is disabled.
      * For nSamplePerBlock > 0, the value is the number of sampled values
      * averaged in each block. 
      */
      int nSamplePerBlock() const;
      
      /**
      * Get number of variables.
      *
      * Call only on processors that have accumulators.
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
      * Call only on processors that have accumulators.
      *
      * \param i integer index of value.
      */
      const Average& accumulator(int i) const;
      
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;  
      
      using ParamComposite::read;
      using Analyzer<D>::interval;
      using Analyzer<D>::isAtInterval;
      using Analyzer<D>::outputFileName;
      

   protected:
   
      using Analyzer<D>::setClassName;
      using Analyzer<D>::readInterval;
      using Analyzer<D>::readOutputFileName;
      
      /**
      * Instantiate Average accumulators and set nSamplePerBlock set nValue.
      *
      * \pre hasAccumulator == false
      * \pre nSamplePerBlock >= 0
      */ 
      void initializeAccumulators(int nValue);
      
      /**
      * Clear internal state of the accumulator.
      *
      * \pre hasAccumulator == true
      */ 
      void clearAccumulators();
      
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
      * Set current value, used by compute function.
      *
      * \param i integer index of variable
      * \param value current value of variable
      */
      void setValue(int i, double value);
      
      /**
      * Compute value of sampled quantity.
      *
      * Call on all processors.
      */
      virtual void compute() = 0;
      
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
      * \pre hasAccumulator() == true
      * \param iStep simulation step counter
      * \param interval analyzer interval, in simulation steps
      */
      void updateAccumulators(long iStep, int interval);

      /**
      * Write results of statistical analysis to files.
      *
      * \pre hasAccumulator == true
      * \param outputFileName base output file name for analyzer
      */
      void outputAccumulators(std::string outputFileName);      
            
      /** 
      * Return reference to parent system.
      */      
      System<D>& system();
      
      // Output file stream
      std::ofstream outputFile_;
      
      // Output filename
      std::string filename_;
      
   private:
      
      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;

      /// Array of current values (only allocated on master processor)
      DArray<double> values_;

      /// Array of value names (only allocated on master processor)
      DArray<std::string> names_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Number of values.
      int nValue_;
 
      /// Does this processor have accumulators ?
      bool hasAccumulators_;


   };
   
   // Inline functions

   /*
   * Get nSamplePerBlock.
   */
   template <int D>
   inline int AverageListAnalyzer<D>::nSamplePerBlock() const
   {  return nSamplePerBlock_; }
   
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
   

}
}
#endif 
