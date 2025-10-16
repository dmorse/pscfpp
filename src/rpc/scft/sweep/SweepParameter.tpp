#ifndef RPC_SWEEP_PARAMETER_TPP
#define RPC_SWEEP_PARAMETER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/system/System.h>
#include <rpc/solvers/Block.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/solvers/MixtureModifier.h>
#include <rpc/solvers/Polymer.h>
#include <rpc/solvers/Solvent.h>
#include <rpc/field/Domain.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/inter/Interaction.h>
#include <pscf/sweep/ParameterModifier.h>
#include <util/global.h>
#include <util/containers/FSArray.h>
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
   SweepParameter<D>::SweepParameter()
    : type_(SweepParameter<D>::Null),
      nId_(0),
      id_(),
      initial_(0.0),
      change_(0.0),
      systemPtr_(nullptr),
      parameterTypesPtr_(nullptr),
      parameterTypeId_(-1)
   {}

   /*
   * Constructor, creates association with system.
   */
   template <int D>
   SweepParameter<D>::SweepParameter(System<D>& system)
    : type_(SweepParameter<D>::Null),
      nId_(0),
      id_(),
      initial_(0.0),
      change_(0.0),
      systemPtr_(&system),
      parameterTypesPtr_(nullptr),
      parameterTypeId_(-1)
   {}

   /*
   * Read type, set nId and allocate id_ array.
   */
   template <int D>
   void SweepParameter<D>::readParamType(std::istream& in)
   {
      std::string buffer;
      in >> buffer;
      std::transform(buffer.begin(), buffer.end(),
                     buffer.begin(), ::tolower);

      if (buffer == "block" || buffer == "block_length") {
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
      } else {
         // Search in parameterTypes array for this sweep parameter
         bool found = false;
         if (parameterTypesPtr_) {
            for (int i = 0; i < parameterTypesPtr_->size(); i++) {
               ParameterType& pType = (*parameterTypesPtr_)[i];
               if (buffer == pType.name) {
                  type_ = Special;
                  nId_ = pType.nId;
                  parameterTypeId_ = i;
                  found = true;
                  break;
               }
            }
         }
         if (!found) {
            std::string msg;
            msg = "Invalid SweepParameter::ParamType value: " + buffer;
            UTIL_THROW(msg.c_str());
         }
      }

      if (id_.isAllocated()) id_.deallocate();
      id_.allocate(nId_);

   }

   /*
   * Write type enum value
   */
   template <int D>
   void SweepParameter<D>::writeParamType(std::ostream& out) const
   {  out << type(); }

   /*
   * Get the ParameterType object for a specialized sweep parameter
   */
   template <int D>
   ParameterType& SweepParameter<D>::parameterType() const
   {  
      UTIL_CHECK(parameterTypesPtr_);
      UTIL_CHECK(isSpecialized());
      return (*parameterTypesPtr_)[parameterTypeId_]; 
   }

   /*
   * Get initial (current) values of swept parameters from parent system.
   */
   template <int D>
   void SweepParameter<D>::getInitial()
   {  initial_ = get_(); }

   /*
   * Set new values of swept parameters in the parent system.
   */
   template <int D>
   void SweepParameter<D>::update(double newVal)
   {  set_(newVal); }

   /*
   * Get string representation of type enum value.
   */
   template <int D>
   std::string SweepParameter<D>::type() const
   {
      if (type_ == Block) {
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
      } else if (type_ == Special) {
         return parameterType().name;
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   double SweepParameter<D>::get_()
   {
      UTIL_CHECK(systemPtr_);
      if (type_ == Chi) {
         return systemPtr_->interaction().chi(id(0), id(1));
      } else if (type_ == Kuhn) {
         return systemPtr_->mixture().monomer(id(0)).kuhn();
      } else if (type_ == Phi_Polymer) {
         return systemPtr_->mixture().polymer(id(0)).phi();
      } else if (type_ == Mu_Polymer) {
         return systemPtr_->mixture().polymer(id(0)).mu();
      } else if (type_ == Block) {
         return systemPtr_->mixture().polymer(id(0)).block(id(1)).length();
      } else if (type_ == Phi_Solvent) {
         return systemPtr_->mixture().solvent(id(0)).phi();
      } else if (type_ == Mu_Solvent) {
         return systemPtr_->mixture().solvent(id(0)).mu();
      } else if (type_ == Solvent) {
         return systemPtr_->mixture().solvent(id(0)).size();
      } else if (type_ == Cell_Param) {
         return systemPtr_->domain().unitCell().parameter(id(0));
      } else if (type_ == Special) {
         ParameterModifier* modifier = parameterType().modifierPtr;
         std::string name = parameterType().name;
         return modifier->getParameter(name,id_);
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   void SweepParameter<D>::set_(double newVal)
   {
      UTIL_CHECK(systemPtr_);
      if (type_ == Chi) {
         systemPtr_->interaction().setChi(id(0), id(1), newVal);
      } else if (type_ == Kuhn) {
         //systemPtr_->mixture().setKuhn(id(0), newVal);
         systemPtr_->mixtureModifier().setKuhn(id(0), newVal);
      } else if (type_ == Phi_Polymer) {
         //systemPtr_->mixture().polymer(id(0)).setPhi(newVal);
         systemPtr_->mixtureModifier().setPhiPolymer(id(0), newVal);
      } else if (type_ == Mu_Polymer) {
         //systemPtr_->mixture().polymer(id(0)).setMu(newVal);
         systemPtr_->mixtureModifier().setMuPolymer(id(0), newVal);
      } else if (type_ == Block) {
         //systemPtr_->mixture().polymer(id(0)).block(id(1)).setLength(newVal);
         systemPtr_->mixtureModifier().setBlockLength(id(0), id(1), newVal);
      } else if (type_ == Phi_Solvent) {
         //systemPtr_->mixture().solvent(id(0)).setPhi(newVal);
         systemPtr_->mixtureModifier().setPhiSolvent(id(0), newVal);
      } else if (type_ == Mu_Solvent) {
         //systemPtr_->mixture().solvent(id(0)).setMu(newVal);
         systemPtr_->mixtureModifier().setMuSolvent(id(0), newVal);
      } else if (type_ == Solvent) {
         //systemPtr_->mixture().solvent(id(0)).setSize(newVal);
         systemPtr_->mixtureModifier().setSolventSize(id(0), newVal);
      } else if (type_ == Cell_Param) {
         FSArray<double,6> params 
                             = systemPtr_->domain().unitCell().parameters();
         params[id(0)] = newVal;
         systemPtr_->setUnitCell(params);
      } else if (type_ == Special) {
         ParameterModifier* modifier = parameterType().modifierPtr;
         std::string name = parameterType().name;
         return modifier->setParameter(name,id_,newVal);
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   template <class Archive>
   void SweepParameter<D>::serialize(Archive ar, const unsigned int version)
   {
      serializeEnum(ar, type_, version);
      ar & nId_;
      for (int i = 0; i < nId_; ++i) {
         ar & id_[i];
      }
      ar & initial_;
      ar & change_;
   }

   // Definitions of operators, with no explicit instantiations.

   /*
   * Inserter for reading a SweepParameter from an istream.
   */
   template <int D>
   std::istream& operator >> (std::istream& in,
                              SweepParameter<D>& param)
   {
      // Read the parameter type identifier string
      param.readParamType(in);

      // Read the identifiers associated with this parameter type.
      for (int i = 0; i < param.nId_; ++i) {
         in >> param.id_[i];
      }
      // Read in the range in the parameter to sweep over
      in >> param.change_;

      return in;
   }

   /*
   * Extractor for writing a SweepParameter to ostream.
   */
   template <int D>
   std::ostream& operator << (std::ostream& out,
                              SweepParameter<D> const & param)
   {
      param.writeParamType(out);
      out << "  ";
      for (int i = 0; i < param.nId_; ++i) {
         out << param.id(i);
         out << " ";
      }
      out << param.change_;

      return out;
   }

}
}

#endif
