#ifndef RPC_TYPES_H
#define RPC_TYPES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>

// Forward declarations
namespace Util {
   template <typename T> class DArray;
   template <typename T> class DRArray;
   class Random;
}
namespace Pscf {
   class CpuVecRandom;
   namespace Prdc {
      class Environment;
      namespace Cpu {
         template <int D> class RField;
         template <int D> class RFieldDft;
         template <int D> class FFT;
         template <int D> class RFieldComparison;
         template <int D> class RFieldDftComparison;
         template <int D> class WaveList;
      }
   }
   namespace Rpc {
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
      template <int D> class EnvironmentFactory;
      //template <int D> class ScftThermo;
      //template <int D> class Iterator;
      //template <int D> class IteratorFactory;
      //template <int D> class Sweep;
      //template <int D> class SweepParameter;
      template <int D> class SweepFactory;
      //template <int D> class BasisFieldState;
      //template <int D> class Simulator;
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
namespace Rpc {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * Aliases for types used in the Rpc program-level namespace.
   *
   * \ingroup Rpc_System_Module
   */
   template <int D>
   class Types
   {

   public:

      using VecRandom = CpuVecRandom;

      template <typename Data> using HostArray = Util::DArray<Data>;
      using Vector = Util::DRArray<double>;

      using Real = double;
      using Complex = fftw_complex;

      using RField = Prdc::Cpu::RField<D>;
      using RFieldDft = Prdc::Cpu::RFieldDft<D>;
      using FFT = Prdc::Cpu::FFT<D>;
      using RFieldComparison = Prdc::Cpu::RFieldComparison<D>;
      using RFieldDftComparison = Prdc::Cpu::RFieldDftComparison<D>;
      using WaveList = Prdc::Cpu::WaveList<D>;

      //using System = Rp::System<D, Rpc::Types<D> >;
      using SystemConstRef = Rpc::SystemConstRef<D>;

      //using Mixture = Rp::Mixture<D, Rpc::Types<D> >;
      //using MixtureModifier = Rp::MixtureModifier<D, Rpc::Types<D> >;
      //using Polymer = Rpc::Polymer<D>;
      //using Solvent = Rpc::Solvent<D>;
      //using Block = Rpc::Block<D>;
      //using Propagator = Propagator<D, Rpc::Types<D> >;

      //using Domain = Rpc::Domain<D>;
      //using FieldIo = Rpc::FieldIo<D>;
      //using WFields = Rpc::WFields<D>;
      //using CFields = Rpc::CFields<D>;
      //using Mask = Rpc::Mask<D>;

      using Environment = Prdc::Environment;
      using EnvironmentFactory = Rpc::EnvironmentFactory<D>;

      //using ScftThermo = Rp::ScftThermo<D, Rpc::Types<D> >;
      //using Iterator = Rp::Iterator<D, Rpc::Types<D> >;
      //using IteratorFactory = Rp::IteratorFactory<D, Rpc::Types<D> >;
      //using Sweep = Rp::Sweep<D, Rpc::Types<D> >;
      //using SweepParameter = Rp::SweepParameter<D, Rpc::Types<D> >;
      //using BasisFieldState = Rp::BasisFieldState<D, Rpc::Types<D> >;
      using SweepFactory = Rpc::SweepFactory<D>;

      //using Simulator = Rp::Simulator<D, Rpc::Types<D> >;
      using SimulatorFactory = Rpc::SimulatorFactory<D>;
      using SimState = Rpc::SimState<D>;

      using Compressor = Rpc::Compressor<D>;
      using CompressorFactory = Rpc::CompressorFactory<D>;
      using IntraCorrelation = Rpc::IntraCorrelation<D>;

      using Perturbation = Rpc::Perturbation<D>;
      using PerturbationFactory = Rpc::PerturbationFactory<D>;

      using Ramp = Rpc::Ramp<D>;
      using RampParameter = Rpc::RampParameter<D>;
      using RampFactory = Rpc::RampFactory<D>;

      using TrajectoryReader = Rpc::TrajectoryReader<D>;
      using TrajectoryReaderFactory = Rpc::TrajectoryReaderFactory<D>;

      using Analyzer = Rpc::Analyzer<D>;
      using AnalyzerFactory = Rpc::AnalyzerFactory<D>;
      using AnalyzerManager = Rpc::AnalyzerManager<D>;
      using AverageAnalyzer = Rpc::AverageAnalyzer<D>;
      using AverageListAnalyzer = Rpc::AverageListAnalyzer<D>;

      using BdSimulator = Rpc::BdSimulator<D>;
      using BdStep = Rpc::BdStep<D>;
      using BdStepFactory = Rpc::BdStepFactory<D>;

      using McSimulator = Rpc::McSimulator<D>;
      using McMove = Rpc::McMove<D>;
      using McMoveFactory = Rpc::McMoveFactory<D>;
      using McMoveManager = Rpc::McMoveManager<D>;

      // Static member functions

      /**
      * Initialize backend.
      */
      static void init();

      /**
      * Set the number of threads.
      *
      * \param nThread  number of threads
      */
      static void setThreadCount(int nThread);

      /**
      * Link vector and scalar random number generators, if needed.
      *
      * The implementation for CPU code associates the two generators.
      * In CPU code, the vector RNG simply calls the scalar RNG.
      *
      * \param vr vector RNG
      * \param sr scalar RNG
      */
      static
      void linkVecRandom(VecRandom& vr, Random & sr);

      /**
      * Set the vector random number generator seed, if needed.
      *
      * The implementation for CPU code does nothing, because the vector
      * RNG is simply a facade for the scalar RNG.
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

} // namespace Rpc
} // namespace Pscf
#endif
