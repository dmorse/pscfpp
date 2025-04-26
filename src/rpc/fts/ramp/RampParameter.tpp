#ifndef RPC_RAMP_PARAMETER_TPP
#define RPC_RAMP_PARAMETER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/perturbation/Perturbation.h>
#include <rpc/solvers/Block.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/Polymer.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/inter/Interaction.h>
#include <util/containers/FSArray.h>
#include <util/global.h>
#include <algorithm>
#include <iomanip>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Default constructor.
   */
   template <int D>
   RampParameter<D>::RampParameter()
    : type_(RampParameter<D>::Null),
      nId_(0),
      id_(),
      initial_(0.0),
      change_(0.0),
      simulatorPtr_(0),
      systemPtr_(0)
   {}

   /*
   * Constructor, creates association with simulator and system.
   */
   template <int D>
   RampParameter<D>::RampParameter(Simulator<D>& simulator)
    : type_(RampParameter<D>::Null),
      nId_(0),
      id_(),
      initial_(0.0),
      change_(0.0),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system()))
   {}


   /*
   * Set the simulator and system associated with this object.
   */
   template <int D>
   void RampParameter<D>::setSimulator(Simulator<D>& simulator)
   {
      simulatorPtr_ = &simulator;
      systemPtr_ = &(simulator.system());
   }

   /*
   * Read type, set nId and allocate id_ array.
   */
   template <int D>
   void RampParameter<D>::readParamType(std::istream& in)
   {
      std::string buffer;
      in >> buffer;
      std::transform(buffer.begin(), buffer.end(),
                     buffer.begin(), ::tolower);

      if (buffer == "block" || buffer == "block_length") {
         UTIL_CHECK(PolymerModel::isThread());
         type_ = Block;
         nId_ = 2; // polymer and block identifiers
      } else if (buffer == "chi") {
         type_ = Chi;
         nId_ = 2; // two monomer type identifiers
      } else if (buffer == "kuhn") {
         type_ = Kuhn;
         nId_ = 1; // monomer type identifier
      } else if (buffer == "phi_polymer") {
         type_ = Phi_Polymer;
         nId_ = 1; // species identifier.
      } else if (buffer == "phi_solvent") {
         type_ = Phi_Solvent;
         nId_ = 1; // species identifier.
      } else if (buffer == "mu_polymer") {
         type_ = Mu_Polymer;
         nId_ = 1; // species identifier.
      } else if (buffer == "mu_solvent") {
         type_ = Mu_Solvent;
         nId_ = 1; // species identifier.
      } else if (buffer == "solvent" || buffer == "solvent_size") {
         type_ = Solvent;
         nId_ = 1; // species identifier.
      } else if (buffer == "cell_param") {
         type_ = Cell_Param;
         nId_ = 1; // lattice parameter identifier.
      } else if (buffer == "lambda_pert") {
         type_ = Lambda_Pert;
         nId_ = 0; // No associated index
      } else if (buffer == "v_monomer") {
         type_ = Vmonomer;
         nId_ = 0; // No associated index
      } else {
         UTIL_THROW("Invalid RampParameter::ParamType value");
      }

      if (id_.isAllocated()) id_.deallocate();
      if (nId_ > 0) {
         id_.allocate(nId_);
      }

   }

   /*
   * Write type enum value
   */
   template <int D>
   void RampParameter<D>::writeParamType(std::ostream& out) const
   {  out << type(); }

   /*
   * Get initial (current) values of swept parameters from parent system.
   */
   template <int D>
   void RampParameter<D>::getInitial()
   {  initial_ = get_(); }

   /*
   * Set new values of swept parameters in the parent system.
   */
   template <int D>
   void RampParameter<D>::update(double newVal)
   {  set_(newVal); }

   /*
   * Get string representation of type enum value.
   */
   template <int D>
   std::string RampParameter<D>::type() const
   {
      if (type_ == Block) {
         UTIL_CHECK(PolymerModel::isThread());
         return "block";
      } else if (type_ == Chi) {
         return "chi";
      } else if (type_ == Kuhn) {
         return "kuhn";
      } else if (type_ == Phi_Polymer) {
         return "phi_polymer";
      } else if (type_ == Phi_Solvent) {
         return "phi_solvent";
      } else if (type_ == Mu_Polymer) {
         return "mu_polymer";
      } else if (type_ == Mu_Solvent) {
         return "mu_solvent";
      } else if (type_ == Solvent) {
         return "solvent_size";
      } else if (type_ == Cell_Param) {
         return "cell_param";
      } else if (type_ == Lambda_Pert) {
         return "lambda_pert";
      } else if (type_ == Vmonomer) {
         return "vMonomer";
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   double RampParameter<D>::get_()
   {
      if (type_ == Block) {
         UTIL_CHECK(PolymerModel::isThread());
         return systemPtr_->mixture().polymer(id(0)).block(id(1)).length();
      } else if (type_ == Chi) {
         return systemPtr_->interaction().chi(id(0), id(1));
      } else if (type_ == Kuhn) {
         return systemPtr_->mixture().monomer(id(0)).kuhn();
      } else if (type_ == Phi_Polymer) {
         return systemPtr_->mixture().polymer(id(0)).phi();
      } else if (type_ == Phi_Solvent) {
         return systemPtr_->mixture().solvent(id(0)).phi();
      } else if (type_ == Mu_Polymer) {
         return systemPtr_->mixture().polymer(id(0)).mu();
      } else if (type_ == Mu_Solvent) {
         return systemPtr_->mixture().solvent(id(0)).mu();
      } else if (type_ == Solvent) {
         return systemPtr_->mixture().solvent(id(0)).size();
      } else if (type_ == Cell_Param) {
         return systemPtr_->domain().unitCell().parameter(id(0));
      } else if (type_ == Lambda_Pert) {
         UTIL_CHECK(simulatorPtr_->hasPerturbation());
         return simulatorPtr_->perturbation().lambda();
      } else if (type_ == Vmonomer) {
         return systemPtr_->mixture().vMonomer();
      }else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   void RampParameter<D>::set_(double newVal)
   {
      if (type_ == Block) {
         UTIL_CHECK(PolymerModel::isThread());
         systemPtr_->mixture().polymer(id(0)).block(id(1)).setLength(newVal);
      } else if (type_ == Chi) {
         systemPtr_->interaction().setChi(id(0), id(1), newVal);
      } else if (type_ == Kuhn) {
         systemPtr_->mixture().setKuhn(id(0), newVal);
      } else if (type_ == Phi_Polymer) {
         systemPtr_->mixture().polymer(id(0)).setPhi(newVal);
      } else if (type_ == Phi_Solvent) {
         systemPtr_->mixture().solvent(id(0)).setPhi(newVal);
      } else if (type_ == Mu_Polymer) {
         systemPtr_->mixture().polymer(id(0)).setMu(newVal);
      } else if (type_ == Mu_Solvent) {
         systemPtr_->mixture().solvent(id(0)).setMu(newVal);
      } else if (type_ == Solvent) {
         systemPtr_->mixture().solvent(id(0)).setSize(newVal);
      } else if (type_ == Cell_Param) {
         FSArray<double,6> params 
                            = systemPtr_->domain().unitCell().parameters();
         UTIL_CHECK(id(0) < params.size());
         params[id(0)] = newVal;
         systemPtr_->setUnitCell(params);
      } else if (type_ == Lambda_Pert) {
         UTIL_CHECK(simulatorPtr_->hasPerturbation());
         return simulatorPtr_->perturbation().setLambda(newVal);
      } else if (type_ == Vmonomer) {
         systemPtr_->mixture().setVmonomer(newVal);
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   template <class Archive>
   void RampParameter<D>::serialize(Archive ar, const unsigned int version)
   {
      serializeEnum(ar, type_, version);
      ar & nId_;
      if (nId_ > 0) {
         for (int i = 0; i < nId_; ++i) {
            ar & id_[i];
         }
      }
      ar & initial_;
      ar & change_;
   }

   // Definitions of operators, with no explicit instantiations.

   /*
   * Inserter for reading a RampParameter from an istream.
   */
   template <int D>
   std::istream& operator >> (std::istream& in,
                              RampParameter<D>& param)
   {
      // Read the parameter type identifier string
      param.readParamType(in);

      // Read the identifiers associated with this parameter type.
      if (param.nId_ > 0) {
         for (int i = 0; i < param.nId_; ++i) {
            in >> param.id_[i];
         }
      }

      // Read in the range in the parameter to sweep over
      in >> param.change_;

      return in;
   }

   /*
   * Extractor for writing a RampParameter to ostream.
   */
   template <int D>
   std::ostream& operator << (std::ostream& out,
                              RampParameter<D> const & param)
   {
      param.writeParamType(out);
      out << "  ";
      if (param.nId_ > 0) {
         for (int i = 0; i < param.nId_; ++i) {
            out << param.id(i);
            out << " ";
         }
      }
      out << param.change_;

      return out;
   }

}
}
#endif
