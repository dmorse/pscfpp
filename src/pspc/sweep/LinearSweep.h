#ifndef PSPC_LINEAR_SWEEP_H
#define PSPC_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"      // base class
#include <util/global.h>
#include <util/containers/DArray.h>

namespace Pscf {
    namespace Pspc {

        template <int D>
        class System;

        using namespace Util;

        /**
         * Base class for a sweep in parameter space where parameters change
         * linearly with the sweep variable. 
         * \ingroup Pspc_Sweep_Module
        */
        template <int D>
        class LinearSweep : public Sweep<D>
        {
        public:

            /** 
             * Constructor.
             * \param system parent System object
             */
            LinearSweep(System<D>& system);

            /**
             * Read parameters from param file.
             */
            void readParameters(std::istream& in);

            /**
             * Setup operation at the beginning of a sweep. Gets initial values of individual parameters.
             */
            void setup();

            /**
             * Set the state before an iteration. Called with each new iteration in SweepTempl::sweep()
             *
             * \param s path length coordinate, in [0,1]
            */    
            void setParameters(double s);

            /**
             * Output data to a running summary.
             *
             * \param out  output file, open for writing
             */
            void outputSummary(std::ostream& out);

            /**
             * Class for storing individual sweep parameters.
             */
            class Parameters;
            
        };

    }
}

#endif