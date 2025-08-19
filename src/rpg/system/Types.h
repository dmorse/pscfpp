#ifndef RPG_TYPES_H
#define RPG_TYPES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declarations
namespace Pscf {
   class Interaction;
   namespace Prdc {
      class Environment;
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
      template <int D> class Mixture;
      template <int D> class MixtureModifier;
      template <int D> class Polymer;
      template <int D> class Solvent;
      template <int D> class Block;
      template <int D> class Propagator;
      template <int D> class Domain;
      template <int D> class FieldIo;
      template <int D> class WFields;
      template <int D> class CFields;
      template <int D> class Mask;
      template <int D> class ScftThermo;
      template <int D> class EnvironmentFactory;
      template <int D> class Iterator;
      template <int D> class IteratorFactory;
      template <int D> class Sweep;
      template <int D> class SweepFactory;
      template <int D> class Simulator;
      template <int D> class SimulatorFactory;
   }
}

namespace Pscf {
namespace Rpg {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * List of aliases for types used the in Rpg namespace.
   *
   * \ingroup Rpg_System_Module
   */
   template <int D>
   class Types
   {

   public:

      using System = Rpg::System<D>;

      using Mixture = Rpg::Mixture<D>;
      using MixtureModifier = Rpg::MixtureModifier<D>;
      using Polymer = Rpg::Polymer<D>;
      using Solvent = Rpg::Solvent<D>;
      using Block = Rpg::Block<D>;
      using Propagator = Rpg::Propagator<D>;

      using Interaction = Pscf::Interaction;
      using Domain = Rpg::Domain<D>;
      using FieldIo = Rpg::FieldIo<D>;
      using ScftThermo = Rpg::ScftThermo<D>;

      using WFields = Rpg::WFields<D>;
      using CFields = Rpg::CFields<D>;
      using Mask = Rpg::Mask<D>;

      using Environment = Prdc::Environment;
      using EnvironmentFactory = Rpg::EnvironmentFactory<D>;
      using Iterator = Rpg::Iterator<D>;
      using IteratorFactory = Rpg::IteratorFactory<D>;
      using Sweep = Rpg::Sweep<D>;
      using SweepFactory = Rpg::SweepFactory<D>;
      using Simulator = Rpg::Simulator<D>;
      using SimulatorFactory = Rpg::SimulatorFactory<D>;

      using RField = Prdc::Cuda::RField<D>;
      using RFieldDft = Prdc::Cuda::RFieldDft<D>;
      using FFT = Prdc::Cuda::FFT<D>;
      using RFieldComparison = Prdc::Cuda::RFieldComparison<D>;
      using RFieldDftComparison = Prdc::Cuda::RFieldDftComparison<D>;
      using WaveList = Prdc::Cuda::WaveList<D>;

   };

} // namespace Rpg
} // namespace Pscf
#endif
