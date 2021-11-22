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
#include <iostream>
#include <iomanip>

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
            class Parameter;

        // friends:
            
        };

        /**
        * Definition of the (eventually private) nested Parameter class, nested in LinearSweep.
        */
        template <int D>
        class LinearSweep<D>::Parameter
        {

        public:
            /**
             * @brief Construct a new Parameter object
             */
            Parameter()
            {}

            void setInitial()
            {

            }
            void update(double s)
            {

            }
            /// Enumeration of allowed parameter types.
            enum paramType { Block, Chi, Kuhn, Phi, Mu, Solvent };
        // friends:
            friend std::istream& operator >> (std::istream&, LinearSweep<D>::Parameter&);
            friend std::ostream& operator << (std::ostream&, LinearSweep<D>::Parameter const&);
            // template <class Archive>
            // friend void serialize(Archive&, LinearSweep<D>::Parameter&, const unsigned int);

        private:
            
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
            /// Change in parameter
            double   change_;
        
        };

        // Declarations of operator templates
        // template <int D>
        // std::istream& operator >> (std::istream& in, typename LinearSweep<D>::Parameter::paramType& type);
        // template <int D>
        // std::ostream& operator << (std::ostream& out, typename LinearSweep<D>::Parameter::paramType type);
        // template <int D>
        // std::istream& operator >> (std::istream& in, typename LinearSweep<D>::Parameter& param);
        // template <int D>
        // std::ostream& operator << (std::ostream& out, typename LinearSweep<D>::Parameter param);

        // Definitions of operators, no explicit instantiations. 
        template <int D>
        std::istream& operator >> (std::istream& in, typename LinearSweep<D>::Parameter::paramType& type)
        {

            std::string buffer;
            in >> buffer;
            std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);

            if (buffer == "block") {
                type = LinearSweep<D>::Parameter::Block;
            } else 
            if (buffer == "chi") {
                type = LinearSweep<D>::Parameter::Chi;
            } else 
            if (buffer == "kuhn") {
                type = LinearSweep<D>::Parameter::Kuhn;
            } else 
            if (buffer == "phi") {
                type = LinearSweep<D>::Parameter::Phi;
            } else 
            if (buffer == "mu") {
                type = LinearSweep<D>::Parameter::Mu;
            } else 
            if (buffer == "solvent") {
                type = LinearSweep<D>::Parameter::Solvent;
            } else {
                UTIL_THROW("Invalid LinearSweep<D>::Parameter::paramType value input");
            }
            return in;
        }

        template <int D>
        std::ostream& operator << (std::ostream& out, typename LinearSweep<D>::Parameter::paramType type)
        {
            if (type == LinearSweep<D>::Parameter::Block) {
                out << "block";
            } else 
            if (type == LinearSweep<D>::Parameter::Chi) {
                out << "chi";
            } else 
            if (type == LinearSweep<D>::Parameter::Kuhn) {
                out << "kuhn";
            } else 
            if (type == LinearSweep<D>::Parameter::Phi) {
                out << "phi";
            } else 
            if (type == LinearSweep<D>::Parameter::Mu) {
                out << "mu";
            } else 
            if (type == LinearSweep<D>::Parameter::Solvent) {
                out << "solvent";
            } else {
                UTIL_THROW("This should never happen.");
            }
            return out;
        }

        template <int D>
        std::istream& operator >> (std::istream& in, typename LinearSweep<D>::Parameter& param)
        {
            in >> param.type_;

            return in;
        }

        template <int D>
        std::ostream& operator << (std::ostream& out, typename LinearSweep<D>::Parameter param)
        {
            out << param.type_;

            return out;
        }


        #ifndef PSPC_LINEAR_SWEEP_TPP
        // Suppress implicit instantiation
        extern template class LinearSweep<1>;
        extern template class LinearSweep<2>;
        extern template class LinearSweep<3>;
        #endif
        
    }
}

#endif