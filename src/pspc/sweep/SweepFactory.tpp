#ifndef PSPC_SWEEP_FACTORY_TPP
#define PSPC_SWEEP_FACTORY_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SweepFactory.h"

// subclasses of Sweep
#include "LinearSweep.h"
// #include "NonlinearSweep.h" Eventually!

namespace Pscf{
    namespace Pspc{

        using namespace Util;

        template <int D>
        SweepFactory<D>::SweepFactory(System<D>& system)
        : systemPtr_(&system)
        {}

        /**
         * Return a pointer to a Sweep subclass with name className
         */
        template <int D>
        Sweep<D>* SweepFactory<D>::factory(std::string const & className) const
        {
            Sweep<D> *ptr = 0;

            // First check if name is known by any subfactories
            ptr = trySubfactories(className);
            if (ptr) return ptr;

            // Explicit class names
            if (className == "LinearSweep") {
                ptr = new LinearSweep<D>(*systemPtr_);
            }

            return ptr;
        }

    }
}


#endif