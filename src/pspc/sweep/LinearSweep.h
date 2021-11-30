#ifndef PSPC_LINEAR_SWEEP_H
#define PSPC_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"      // base class
#include <pscf/inter/ChiInteraction.h>
#include <util/global.h>
#include <iostream>
#include <algorithm>
#include <iomanip>

namespace Pscf {
namespace Pspc {

   // Forward declare classes and operators
   template <int D>
   class System;

   template <int D>
   class LinearSweep;

   template <int D>
   class LinearSweepParameter;

   using namespace Util;

   /**
   * Base class for a sweep in parameter space where parameters change
   * linearly with the sweep variable. 
   * 
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
      * 
      * \param in Input stream from param file.
      */
      void readParameters(std::istream& in);

      /**
      * Setup operation at the beginning of a sweep. Gets initial 
      * values of individual parameters.
      */
       void setup();

      /**
      * Set the state before an iteration. Called with each new iteration 
      * in SweepTempl::sweep()
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

   protected:

      using Sweep<D>::system;
      using Sweep<D>::hasSystem;
   
   private:

      /**
      * Number of parameters being swept. 
      */
      int nParameter_; 

      /**
      * Array of LinearSweepParameter objects.
      */
      DArray< LinearSweepParameter<D> > parameters_;

   // friends:
       
   };

   /**
   * Class for storing data about an individual sweep parameter.
   * 
   * \ingroup Pspc_Sweep_Module
   */
   template <int D>
   class LinearSweepParameter
   {

   public:

      /**
      * Construct a new LinearSweepParameter object.
      */
      LinearSweepParameter()
       : id_(),
         systemPtr_(0)
      {}

      /**
      * Construct a LinearSweepParameter object with a pointer to the system.
      */
      LinearSweepParameter(System<D>& system)
       : id_(),
         systemPtr_(&system)
      {}

      /**
      * Read type of parameter being swept, and set number of identifiers.
      * 
      * \param in Input stream from param file. 
      */
      void readParamType(std::istream& in)
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
            UTIL_THROW("I don't know how to implement 'solvent'.");
         } else {
            UTIL_THROW("Invalid LinearSweepParameter::paramType value input");
         }

         id_.allocate(nID_);
      }

      /**
      * Write type of parameter swept.
      * 
      * \param out Output file stream.
      */
      void writeParamType(std::ostream& out) const
      {
         out << type();
      }

      /**
      * Store the pre-sweep value of the corresponding parameter.
      */
      void getInitial()
      {
         initial_ = get_();
      }

      /**
      * Update the corresponding parameter value in the system. 
      * 
      * \param s sweep coordinate, with range [0,1]
      */
      void update(double s)
      {
         set_(initial_+s*change_);
      }
      
      /**
      * Return the current system parameter value. 
      */
      double current()
      {
         return get_();
      }

      /**
      * Return a string describing the parameter type for this object. 
      */
      std::string type() const
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

      /**
      * Accessor for the id_ array.
      * 
      * \param i array index to access
      */
      int id(int i) const
      {
         return id_[i];
      }
      
      /**
      * Return the total change planned for this parameter during sweep.
      */
      double change() const
      {
         return change_;
      }

      /**
      * Set the system associated with this object.
      */
      void setSystem(System<D>& system)
      {
         systemPtr_ = &system;
      }
      
      /**
      * Serialize to or from an archive.
      *
      * \param ar Archive object 
      * \param version archive format version index
      */
      template <class Archive>
      void serialize(Archive ar, const unsigned int version)
      {
         serializeEnum(ar, type_, version);
         ar & nID_;
         for (int i = 0; i < nID_; ++i) {
            ar & id_[i];
         }
         ar & initial_;
         ar & change_;
      }

   private:

      /// Enumeration of allowed parameter types.
      enum paramType { Block, Chi, Kuhn, Phi_Polymer, Phi_Solvent,
                                     Mu_Polymer, Mu_Solvent, Solvent};

      /// Type of parameter associated with an object of this class. 
      paramType type_;
      
      /// Number of identifiers needed for this parameter type. 
      int nID_;

      /// Identifier indices.
      DArray<int> id_;

      /// Initial value of parameter 
      double   initial_;

      /// Change in parameter
      double   change_;
      
      /// Pointer to the overall system. Set in LinearSweepParameter constructor.
      System<D>* systemPtr_;

      // Gets the current system parameter value.
      double get_()
      {
         if (type_ == Block) {
            return systemPtr_->mixture().polymer(id(0)).block(id(1)).length();
         } else 
         if (type_ == Chi) {
            return systemPtr_->interaction().chi(id(0),id(1));
         } else 
         if (type_ == Kuhn) {
            return systemPtr_->mixture().monomer(id(0)).step();
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

      // Sets the system parameter value. 
      void set_(double newVal)
      {
         if (type_ == Block) {
            systemPtr_->mixture().polymer(id(0)).block(id(1)).setLength(newVal);
         } else 
         if (type_ == Chi) {
            systemPtr_->interaction().setChi(id(0),id(1),newVal);
         } else 
         if (type_ == Kuhn) {
            int monomerId;
            double kuhn;
            // Set new value of kuhn length for this monomer.
            systemPtr_->mixture().monomer(id(0)).setStep(newVal);
            // Update the kuhn lengths of all blocks. Potential mismatch otherwise.
            for (int i=0; i < systemPtr_->mixture().nPolymer(); ++i) {
               for (int j=0; j < systemPtr_->mixture().polymer(i).nBlock(); ++j) {
                  monomerId = systemPtr_->mixture().polymer(i).block(j).monomerId();
                  kuhn = systemPtr_->mixture().monomer(monomerId).step();
                  systemPtr_->mixture().polymer(i).block(j).setKuhn(kuhn);
               }
            }
            // Solvent kuhn doesn't need updating. 
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

   // friends:
      template <int U>
      friend std::istream& operator >> (std::istream&, LinearSweepParameter<U>&);

      template <int U>
      friend std::ostream& operator << (std::ostream&, LinearSweepParameter<U> const&);
   
   };

   // Definitions of operators, with no explicit instantiations. 

   /**
   * Inserter operator for reading data into a LinearSweepParameter.
   */
   template <int D>
   std::istream& operator >> (std::istream& in, 
                              LinearSweepParameter<D>& param)
   {
      // read the parameter type.
      param.readParamType(in);  
      // read the identifiers associated with this parameter type. 
      for (int i = 0; i < param.nID_; ++i) {
         in >> param.id_[i];
      }
      // read in the range in the parameter to sweep over.
      in >> param.change_;

      return in;
   }

   /**
   * Extractor operator for outputting data from a LinearSweepParameter.
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, 
                              LinearSweepParameter<D> const & param)
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

   #ifndef PSPC_LINEAR_SWEEP_TPP
   // Suppress implicit instantiation
   extern template class LinearSweep<1>;
   extern template class LinearSweep<2>;
   extern template class LinearSweep<3>;
   #endif

}
}
#endif
