#ifndef PSPC_LINEAR_SWEEP_TPP
#define PSPC_LINEAR_SWEEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/System.h>
#include "LinearSweep.h"

namespace Pscf {
    namespace Pspc {

        using namespace Util;

        template <int D>
        LinearSweep<D>::LinearSweep(System<D>& system)
         : Sweep<D>(system)
        {}

        template <int D>
        void LinearSweep<D>::readParameters(std::istream& in)
        {}

        template <int D>
        void LinearSweep<D>::setup()
        {}

        template <int D>
        void LinearSweep<D>::setParameters(double s)
        {}

        template <int D>
        void LinearSweep<D>::outputSummary(std::ostream& out)
        {}

    }
}

#endif