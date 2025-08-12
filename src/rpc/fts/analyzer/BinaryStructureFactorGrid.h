#ifndef RPC_BINARY_STRUCTURE_FACTOR_GRID_H
#define RPC_BINARY_STRUCTURE_FACTOR_GRID_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                             // base class

#include <util/containers/DArray.h>               // member
#include <util/accumulators/Average.h>            // member
#include <prdc/cpu/RField.h>                      // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <map>                                    // member

#include <string>
#include <iostream>
#include <vector>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * BinaryStructureFactorGrid evaluates AB copolymer structure factors.
   *
   * This class evaluates the structures factors for all wavevectors within
   * a specified region of a Fourier space grid.
   * 
   * A structure factor for a wavevector k for AB diblock defined as an 
   * expectation value
   * \f[
   *     S(k)  = n/(V \chi N)^2 <W_(k)W_(-k)> - 1/(2 \chi N)
   * \f]
   * where, V is system volume, and \f$W_(k)\f$ is a Fourier mode of 
   * fluctuating field 
   *
   * \see \ref rpc_BinaryStructureFactorGrid_page "Manual Page"
   * 
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryStructureFactorGrid : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      */
      BinaryStructureFactorGrid(Simulator<D>& simulator, System<D>& system);

      /**	
      * Destructor.
      */
      ~BinaryStructureFactorGrid(){};

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
      * Compute structure factor
      */
      void computeStructureFactor();
      
      /**
      * Compute average S(k) over k of equal magnitude
      */
      void averageStructureFactor();
   
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
      
      /**
      * Structure factor
      */
      DArray<double> structureFactors_;
   
      /// Number of wavevectors.
      int nWave_;

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
      using ParamComposite::readOptional;
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

      /// Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in wavevector mesh 
      int kSize_;

      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;
      
      /// wField in Fourier mode
      DArray< RFieldDft<D> > wKGrid_;

      /// W_ for diblock copolymer (wm = (wa-wb)/2)
      RField<D> wm_;

      /// W_ in Fourier mode
      RFieldDft<D> wk_;

      /// Bare wavenumber value q = sqrt(kSq) lists
      std::vector<double> qList_;
      
      /// Map key to be qsquare and value to be average structure factor over k of equal magnitude
      std::map<double, double> averageSMap_;

   };
   
   // Get the parent system.
   template <int D>
   inline System<D>& BinaryStructureFactorGrid<D>::system()
   {  return *systemPtr_; }
   
   //Get parent Simulator object.
   template <int D>
   inline Simulator<D>& BinaryStructureFactorGrid<D>::simulator()
   {  return *simulatorPtr_; }

   #ifndef RPC_BINARY_STRUCTURE_FACTOR_GRID_TPP
   // Suppress implicit instantiation
   extern template class BinaryStructureFactorGrid<1>;
   extern template class BinaryStructureFactorGrid<2>;
   extern template class BinaryStructureFactorGrid<3>;
   #endif

}
}
#endif
