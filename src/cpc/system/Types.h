#ifndef CPC_TYPES_H
#define CPC_TYPES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declarations
namespace Pscf {
   class Interaction;
   template <typename WT> class Species;
   template <typename WT> class PolymerSpecies;
   template <typename WT> class SolventSpecies;
   template <typename WT> class Composition;
   namespace Prdc {
      class Environment;
      namespace Cpu {
         //template <int D> class CField;
         //template <int D> class FFT;
         //template <int D> class CFieldComparison;
         //template <int D> class WaveList;
      }
   }
   namespace Cpc {
      template <int D> class System;
      template <int D> class Mixture;
      // template <int D> class MixtureModifier;
      template <int D> class Polymer;
      template <int D> class Solvent;
      template <int D> class Block;
      template <int D> class Propagator;
      template <int D> class Domain;
      template <int D> class FieldIo;
      template <int D> class WFields;
      template <int D> class CFields;
      // template <int D> class Simulator;
      // template <int D> class SimulatorFactory;
   }
}

namespace Pscf {
namespace Cpc {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Prdc;

   /**
   * List of aliases for types used the in Cpc namespace.
   *
   * \ingroup Cpc_System_Module
   */
   template <int D>
   class Types
   {

   public:

      using System = Cpc::System<D>;

      using Mixture = Cpc::Mixture<D>;
      // using MixtureModifier = Cpc::MixtureModifier<D>;
      using Polymer = Cpc::Polymer<D>;
      using Solvent = Cpc::Solvent<D>;
      using Block = Cpc::Block<D>;
      using Propagator = Cpc::Propagator<D>;

      using Species = Pscf::Species<double>;
      using PolymerSpecies = Pscf::PolymerSpecies<double>;
      using SolventSpecies = Pscf::SolventSpecies<double>;
      using Composition = Pscf::Composition<double>;

      using Interaction = Pscf::Interaction;
      using Domain = Cpc::Domain<D>;
      using FieldIo = Cpc::FieldIo<D>;

      using WFields = Cpc::WFields<D>;
      using CFields = Cpc::CFields<D>;

      // using Simulator = Cpc::Simulator<D>;
      // using SimulatorFactory = Cpc::SimulatorFactory<D>;

      //using CField = Prdc::CField<D,CPT>;
      //using FFT = Prdc::FFT<D,CPT>;
      //using CFieldComparison = Prdc::CFieldComparison<D,CPT>;
      //using WaveList = Prdc::WaveList<D,CPT>;

   };

} // namespace Cpc
} // namespace Pscf
#endif
