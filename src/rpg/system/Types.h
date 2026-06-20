#ifndef RPG_TYPES_H
#define RPG_TYPES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>

// Forward declarations
namespace Pscf {
   class CudaVecRandom;
   template <typename T> class DeviceArray;
   template <typename T> class HostDArray;
   namespace Prdc {
      namespace Cuda {
         template <int D> class RField;
         template <int D> class RFieldDft;
         template <int D> class FFT;
         template <int D> class RFieldComparison;
         template <int D> class RFieldDftComparison;
         template <int D> class WaveList;
      }
   }
   namespace Rpg {
      template <int D> class System;
      template <int D> class SystemConstRef;
      template <int D> class Mixture;
      template <int D> class MixtureModifier;
      //template <int D> class Polymer;
      //template <int D> class Solvent;
      template <int D> class Block;
      template <int D> class Propagator;
      //template <int D> class Domain;
      //template <int D> class FieldIo;
      //template <int D> class WFields;
      //template <int D> class CFields;
      //template <int D> class Mask;
      template <int D> class ScftThermo;
      template <int D> class EnvironmentFactory;
      template <int D> class Iterator;
      template <int D> class IteratorFactory;
      template <int D> class Sweep;
      template <int D> class SweepParameter;
      template <int D> class BasisFieldState;
      template <int D> class SweepFactory;
      template <int D> class Simulator;
      template <int D> class SimulatorFactory;
      template <int D> struct SimState;
      template <int D> class Compressor;
      template <int D> class CompressorFactory;
      template <int D> class IntraCorrelation;
      template <int D> class Perturbation;
      template <int D> class PerturbationFactory;
      template <int D> class Ramp;
      template <int D> class RampParameter;
      template <int D> class RampFactory;
      template <int D> class Analyzer;
      template <int D> class AnalyzerFactory;
      template <int D> class AnalyzerManager;
      template <int D> class AverageAnalyzer;
      template <int D> class AverageListAnalyzer;
      template <int D> class TrajectoryReader;
      template <int D> class TrajectoryReaderFactory;
      template <int D> class BdSimulator;
      template <int D> class BdStep;
      template <int D> class BdStepFactory;
      template <int D> class McSimulator;
      template <int D> class McMove;
      template <int D> class McMoveFactory;
      template <int D> class McMoveManager;
   }
}

namespace Pscf {
namespace Rpg {

   // Namespaces that may be used implicitly
   using namespace Prdc;

   /**
   * List of aliases for types used in the Rpg program-level namespace.
   *
   * \ingroup Rpg_System_Module
   */
   template <int D>
   class Types {

   public:

      using VecRandom = Pscf::CudaVecRandom;

      template <typename T> using HostArray = Pscf::HostDArray<T>;
      using Vector = Pscf::DeviceArray<Pscf::cudaReal>;

      using Real = Pscf::cudaReal;
      using Complex = Pscf::cudaComplex;

      using RField = Prdc::Cuda::RField<D>;
      using RFieldDft = Prdc::Cuda::RFieldDft<D>;
      using FFT = Prdc::Cuda::FFT<D>;
      using RFieldComparison = Prdc::Cuda::RFieldComparison<D>;
      using RFieldDftComparison = Prdc::Cuda::RFieldDftComparison<D>;
      using WaveList = Prdc::Cuda::WaveList<D>;

      using System = Rpg::System<D>;
      using SystemConstRef = Rpg::SystemConstRef<D>;

      using Mixture = Rpg::Mixture<D>;
      using MixtureModifier = Rpg::MixtureModifier<D>;
      //using Polymer = Rpg::Polymer<D>;
      //using Solvent = Rpg::Solvent<D>;
      using Block = Rpg::Block<D>;
      using Propagator = Rpg::Propagator<D>;

      //using Domain = Rpg::Domain<D>;
      //using FieldIo = Rpg::FieldIo<D>;
      //using WFields = Rpg::WFields<D>;
      //using CFields = Rpg::CFields<D>;
      //using Mask = Rpg::Mask<D>;
      using EnvironmentFactory = Rpg::EnvironmentFactory<D>;

      using ScftThermo = Rpg::ScftThermo<D>;

      using Iterator = Rpg::Iterator<D>;
      using IteratorFactory = Rpg::IteratorFactory<D>;
      using Sweep = Rpg::Sweep<D>;
      using SweepParameter = Rpg::SweepParameter<D>;
      using BasisFieldState = Rpg::BasisFieldState<D>;
      using SweepFactory = Rpg::SweepFactory<D>;

      using Simulator = Rpg::Simulator<D>;
      using SimulatorFactory = Rpg::SimulatorFactory<D>;
      using SimState = Rpg::SimState<D>;

      using Compressor = Rpg::Compressor<D>;
      using CompressorFactory = Rpg::CompressorFactory<D>;
      using IntraCorrelation = Rpg::IntraCorrelation<D>;

      using Perturbation = Rpg::Perturbation<D>;
      using PerturbationFactory = Rpg::PerturbationFactory<D>;

      using Ramp = Rpg::Ramp<D>;
      using RampParameter = Rpg::RampParameter<D>;
      using RampFactory = Rpg::RampFactory<D>;

      using Analyzer = Rpg::Analyzer<D>;
      using AnalyzerFactory = Rpg::AnalyzerFactory<D>;
      using AnalyzerManager = Rpg::AnalyzerManager<D>;
      using AverageAnalyzer = Rpg::AverageAnalyzer<D>;
      using AverageListAnalyzer = Rpg::AverageListAnalyzer<D>;

      using TrajectoryReader = Rpg::TrajectoryReader<D>;
      using TrajectoryReaderFactory = Rpg::TrajectoryReaderFactory<D>;

      using BdSimulator = Rpg::BdSimulator<D>;
      using BdStep = Rpg::BdStep<D>;
      using BdStepFactory = Rpg::BdStepFactory<D>;

      using McSimulator = Rpg::McSimulator<D>;
      using McMove = Rpg::McMove<D>;
      using McMoveFactory = Rpg::McMoveFactory<D>;
      using McMoveManager = Rpg::McMoveManager<D>;

      // Static member functions

      /**
      * Initialize backend thread array.
      */
      static void init();

      /**
      * Set the number of threads per block.
      *
      * \param nThread  number of threads
      */
      static void setThreadCount(int nThread);

   };

} // namespace Rpg
} // namespace Pscf
#endif
