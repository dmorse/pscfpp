#ifndef RPC_SYSTEM_CPP
#define RPC_SYSTEM_CPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <prdc/system/SystemTmpl.tpp>

#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/scft/ScftThermo.h>
#include <rpc/environment/EnvironmentFactory.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/simulator/SimulatorFactory.h>
#include <rpc/fts/compressor/Compressor.h>
#include <rpc/scft/iterator/Iterator.h>
#include <rpc/scft/iterator/IteratorFactory.h>
#include <rpc/scft/sweep/Sweep.h>
#include <rpc/scft/sweep/SweepFactory.h>

#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/environment/Environment.h>

#include <pscf/inter/Interaction.h>

// Explicit template instantiation
namespace Pscf {
   namespace Prdc {
      template class SystemTmpl< 1, Rpc::Types<1> >;
      template class SystemTmpl< 2, Rpc::Types<2> >;
      template class SystemTmpl< 3, Rpc::Types<3> >;
   }
   namespace Rpc {

      template <int D>
      System<D>::System()
       : SystemTmpl<D, Types<D> >(*this)
      {  ParamComposite::setClassName("System"); }

      template class System<1>;
      template class System<2>;
      template class System<3>;
   } 
} 
#endif
