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
namespace Util {
   class Random;
}
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
      //template <int D> class System;
      template <int D> class SystemConstRef;
      //template <int D> class Mixture;
      //template <int D> class MixtureModifier;
      //template <int D> class Polymer;
      //template <int D> class Solvent;
      //template <int D> class Block;
      //template <int D> class Propagator;
      //template <int D> class Domain;
      //template <int D> class FieldIo;
      //template <int D> class WFields;
      //template <int D> class CFields;
      //template <int D> class Mask;
      //template <int D> class ScftThermo;
      template <int D> class EnvironmentFactory;
      //template <int D> class Iterator;
      //template <int D> class IteratorFactory;
      //template <int D> class Sweep;
      //template <int D> class SweepParameter;
      //template <int D> class BasisFieldState;

      template <int D> class SweepFactory;

      //template <int D> class Simulator;
      //template <int D> struct SimState;
      //template <int D> class SimulatorFactory;

      template <int D> class Compressor;
      template <int D> class CompressorFactory;
      template <int D> class IntraCorrelation;

      template <int D> class Analyzer;
      template <int D> class AnalyzerFactory;
      template <int D> class AnalyzerManager;
      template <int D> class AverageAnalyzer;
      template <int D> class AverageListAnalyzer;

      //template <int D> class BdSimulator;
      //template <int D> class BdStep;
      template <int D> class BdStepFactory;

      //template <int D> class McSimulator;
      //template <int D> class McMove;
      //template <int D> class McMoveManager;
      template <int D> class McMoveFactory;

      template <int D> class TrajectoryReader;
      template <int D> class TrajectoryReaderFactory;

      template <int D> class Perturbation;
      template <int D> class PerturbationFactory;

      template <int D> class Ramp;
      template <int D> class RampParameter;
      template <int D> class RampFactory;
   }
}

namespace Pscf {
namespace Rpg {

   // Namespaces that may be used implicitly
   using namespace Util;
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

      //using System = Rp::System<D, Rpg::Types<D> >;
      using SystemConstRef = Rpg::SystemConstRef<D>;

      //using Mixture = Rp::Mixture<D, Rpg::Types<D> >;
      //using MixtureModifier = Rp::MixtureModifier<D, Rpg::Types<D> >;
      //using Polymer = Rpg::Polymer<D>;
      //using Solvent = Rpg::Solvent<D>;
      //using Block = Rpg::Block<D>;
      //using Propagator = Rpg::Propagator<D>;

      //using Domain = Rpg::Domain<D>;
      //using FieldIo = Rpg::FieldIo<D>;
      //using WFields = Rpg::WFields<D>;
      //using CFields = Rpg::CFields<D>;
      //using Mask = Rpg::Mask<D>;

      using EnvironmentFactory = Rpg::EnvironmentFactory<D>;

      //using ScftThermo = Rp::ScftThermo<D, Rpg::Types<D> >;
      //using Iterator = Rp::Iterator<D, Rpg::Types<D> >;
      //using IteratorFactory = Rp::IteratorFactory<D, Rpg::Types<D> >;
      //using Sweep = Rp::Sweep<D, Rpg::Types<D> >;
      //using SweepParameter = Rp::SweepParameter<D, Rpg::Types<D> >;
      //using BasisFieldState = Rp::BasisFieldState<D, Rpg::Types<D> >;
      using SweepFactory = Rpg::SweepFactory<D>;

      //using Simulator = Rp::Simulator<D, Rpg::Types<D> >;
      //using SimState = Rp::SimState<D, Rpg::Types<D> >;
      //using SimulatorFactory = Rp::SimulatorFactory<D, Rpg::Types<D> >;

      using Compressor = Rpg::Compressor<D>;
      using CompressorFactory = Rpg::CompressorFactory<D>;
      using IntraCorrelation = Rpg::IntraCorrelation<D>;

      //using BdSimulator = Rp::BdSimulator<D, Rpg::Types<D> >;
      //using BdStep = Rp::BdStep<D, Rpg::Types<D> >;
      using BdStepFactory = Rpg::BdStepFactory<D>;

      //using McSimulator = Rp::McSimulator<D, Rpg::Types<D> >;
      //using McMove = Rp::McMove<D, Rpg::Types<D> >;
      //using McMoveManager = Rp::McMoveManager<D, Rpg::Types<D> >;
      using McMoveFactory = Rpg::McMoveFactory<D>;

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

      /**
      * Link vector and scalar random number generators, if needed.
      *
      * The Cuda GPU implementation does nothing, because the host and
      * device RNGs are independent.
      *
      * \param vr  vector RNG
      * \param sr  scalar RNG
      */
      static
      void linkVecRandom(VecRandom& vr, Random & sr);

      /**
      * Set the vector random number generator seed, if needed.
      *
      * The Cuda GPU implementation sets a seed for the device RNG.
      *
      * \param vr  vector RNG to be initialized
      * \param seed  random seed value
      */
      static
      void seedVecRandom(VecRandom& vr, long seed);

   };

   // Explicit instantiation declarations
   extern template class Types<1>;
   extern template class Types<2>;
   extern template class Types<3>;

} // namespace Rpg
} // namespace Pscf
#endif
