#ifndef PSPC_SWEEP_PARAMETER_TPP
#define PSPC_SWEEP_PARAMETER_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/solvers/Block.h>
#include <pspc/solvers/Mixture.h>
#include <pspc/solvers/Polymer.h>
#include <pspc/System.h>
#include <pscf/inter/ChiInteraction.h>
#include <util/global.h>
#include <algorithm>
#include <iomanip>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   SweepParameter<D>::SweepParameter()
    : type_(SweepParameter<D>::Null),
      nID_(0),
      id_(),
      initial_(0.0),
      change_(0.0),
      systemPtr_(0)
   {}

   /*
   * Constructor, creates association with system.
   */
   template <int D>
   SweepParameter<D>::SweepParameter(System<D>& system)
    : type_(SweepParameter<D>::Null),
      nID_(0),
      id_(),
      initial_(0.0),
      change_(0.0),
      systemPtr_(&system)
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
         nID_ = 2; // polymer and block identifiers
      } else 
      if (buffer == "chi") {
         type_ = Chi;
         nID_ = 2; // two monomer type identifiers
      } else 
      if (buffer == "kuhn") {
         type_ = Kuhn;
         nID_ = 1; // monomer type identifier
      } else 
      if (buffer == "phi_polymer") {
         type_ = Phi_Polymer;
         nID_ = 1; //species identifier.
      } else 
      if (buffer == "phi_solvent") {
         type_ = Phi_Solvent;
         nID_ = 1; //species identifier.
      } else 
      if (buffer == "mu_polymer") {
         type_ = Mu_Polymer;
         nID_ = 1; //species identifier.
      } else 
      if (buffer == "mu_solvent") {
         type_ = Mu_Solvent;
         nID_ = 1; //species identifier.
      } else 
      if (buffer == "solvent" || buffer == "solvent_size") {
         type_ = Solvent;
         nID_ = 1; //species identifier.
      } else {
         UTIL_THROW("Invalid SweepParameter::ParamType value");
      }

      if (id_.isAllocated()) id_.deallocate();
      id_.allocate(nID_);

   }

   /*
   * Write type enum value
   */
   template <int D>
   void SweepParameter<D>::writeParamType(std::ostream& out) const
   {
      out << type();
   }

   /*
   * Get the current value from the parent system.
   */
   template <int D>
   void SweepParameter<D>::getInitial()
   {
      initial_ = get_();
   }

   /*
   * Set a new value in the parent system.
   */
   template <int D>
   void SweepParameter<D>::update(double newVal)
   {
      set_(newVal);
   }

   /*
   * Get string representation of type enum value.
   */
   template <int D>
   std::string SweepParameter<D>::type() const
   {
      if (type_ == Block) {
         return "block_length";
      } else 
      if (type_ == Chi) {
         return "chi";
      } else 
      if (type_ == Kuhn) {
         return "kuhn";
      } else 
      if (type_ == Phi_Polymer) {
         return "phi_polymer";
      } else 
      if (type_ == Phi_Solvent) {
         return "phi_solvent";
      } else  
      if (type_ == Mu_Polymer) {
         return "mu_polymer";
      } else
      if (type_ == Mu_Solvent) {
         return "mu_solvent";
      } else  
      if (type_ == Solvent) {
         return "solvent_size";
      } else {
         UTIL_THROW("This should never happen.");
         }
   }

   // template <int D>
   // SweepParameter<D>::

   template <int D>
   double SweepParameter<D>::get_()
   {
      if (type_ == Block) {
         return systemPtr_->mixture().polymer(id(0)).block(id(1)).length();
      } else 
      if (type_ == Chi) {
         return systemPtr_->interaction().chi(id(0),id(1));
      } else 
      if (type_ == Kuhn) {
         return systemPtr_->mixture().monomer(id(0)).kuhn();
      } else 
      if (type_ == Phi_Polymer) {
         return systemPtr_->mixture().polymer(id(0)).phi();
      } else
      if (type_ == Phi_Solvent) {
         return systemPtr_->mixture().solvent(id(0)).phi();
      } else 
      if (type_ == Mu_Polymer) {
         return systemPtr_->mixture().polymer(id(0)).mu();
      } else
      if (type_ == Mu_Solvent) {
         return systemPtr_->mixture().solvent(id(0)).mu();
      } else 
      if (type_ == Solvent) {
         return systemPtr_->mixture().solvent(id(0)).size();
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   void SweepParameter<D>::set_(double newVal)
   {
      if (type_ == Block) {
         systemPtr_->mixture().polymer(id(0)).block(id(1)).setLength(newVal);
      } else 
      if (type_ == Chi) {
         systemPtr_->interaction().setChi(id(0), id(1), newVal);
      } else 
      if (type_ == Kuhn) {
         Pscf::Pspc::Mixture<D>& mixture = systemPtr_->mixture();

         // Set new kuhn length for this monomer
         mixture.monomer(id(0)).setKuhn(newVal);

         // Update kuhn length for all blocks of this monomer type
         Pscf::Pspc::Block<D>* blockPtr;
         for (int i=0; i < mixture.nPolymer(); ++i) {
            for (int j=0; j < mixture.polymer(i).nBlock(); ++j) {
               blockPtr = &(mixture.polymer(i).block(j));
               if (id(0) == blockPtr->monomerId()) {
                  blockPtr->setKuhn(newVal);
               }
            }
         }

      } else
      if (type_ == Phi_Polymer) {
         systemPtr_->mixture().polymer(id(0)).setPhi(newVal);
      } else
      if (type_ == Phi_Solvent) {
         systemPtr_->mixture().solvent(id(0)).setPhi(newVal);
      } else 
      if (type_ == Mu_Polymer) {
         systemPtr_->mixture().polymer(id(0)).setMu(newVal);
      } else
      if (type_ == Mu_Solvent) {
         systemPtr_->mixture().solvent(id(0)).setMu(newVal);
      } else 
      if (type_ == Solvent) {
         systemPtr_->mixture().solvent(id(0)).setSize(newVal);
      } else {
         UTIL_THROW("This should never happen.");
      }
   }

   template <int D>
   template <class Archive>
   void SweepParameter<D>::serialize(Archive ar, const unsigned int version)
   {
      serializeEnum(ar, type_, version);
      ar & nID_;
      for (int i = 0; i < nID_; ++i) {
         ar & id_[i];
      }
      ar & initial_;
      ar & change_;
   }

   // Definitions of operators, with no explicit instantiations. 

   /**
   * Inserter for reading a SweepParameter from an istream.
   *
   * \param in  input stream
   * \param param  SweepParameter<D> object to read
   */
   template <int D>
   std::istream& operator >> (std::istream& in, 
                              SweepParameter<D>& param)
   {
      // Read the parameter type.
      param.readParamType(in);  
      // Read the identifiers associated with this parameter type. 
      for (int i = 0; i < param.nID_; ++i) {
         in >> param.id_[i];
      }
      // Read in the range in the parameter to sweep over
      in >> param.change_;

      return in;
   }

   /**
   * Extractor for writing a SweepParameter to ostream.
   *
   * \param out  output stream
   * \param param  SweepParameter<D> object to write
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, 
                              SweepParameter<D> const & param)
   {
      param.writeParamType(out);
      out << "  ";
      for (int i = 0; i < param.nID_; ++i) {
         out << param.id(i);
         out << " ";
      }
      out << param.change_;

      return out;
   }

}
}

#endif
