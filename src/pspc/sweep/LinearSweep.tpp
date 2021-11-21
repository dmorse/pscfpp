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

        /**
        * Definition of the SweepParam class nested within Sweep.
        */
        template <int D>
        class LinearSweep<D>::Parameters
        {
            /**
             * @brief Construct a new Parameters object
             */
            Parameters()
            {}

        public:
            void setInitial()
            {
                
            }
            void update(double s)
            {

            }

        private:
            /// Enumeration of allowed parameter types.
            enum paramType { Block, Chi, Kuhn, Phi, Mu, Solvent };
            /// Type of parameter associated with an object of this class. 
            paramType type_;
            /// Number of species associated with parameter type. For example, Kuhn length is of 1 species, Chi is between two.
            int nID_;
            /// List of identifiers of species associated with the parameter type, of length nID.
            DArray<int> id_;
            /// Initial value of parameter 
            double   initial_;
            /// Final value of parameter
            double   final_;
            /// Change in parameter (??? do we want to do this, or the initial, step size, and number of steps as in FORTRAN pscf)
            double   change_;

        // friends:
            //template <int D>
            //friend std::istream& operator >> (std::istream&, Sweep<D>::SweepParam& param)

        };        
        
    }
}

#endif