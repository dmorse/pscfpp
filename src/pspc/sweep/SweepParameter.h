#ifndef PSPC_SWEEP_PARAMETER_H
#define PSPC_SWEEP_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/solvers/Block.h>
#include <pspc/solvers/Mixture.h>
#include <pspc/solvers/Polymer.h>
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

   /**
   * Class for storing data about an individual sweep parameter.
   * 
   * \ingroup Pspc_Sweep_Module
   */
   template <int D>
   class SweepParameter
   {

   public:

      /**
      * Default constructor.
      */
      SweepParameter();

      /**
      * Constructor that stores a pointer to parent system.
      *
      * \param system  parent system
      */
      SweepParameter(System<D>& system);

      /**
      * Read type of parameter being swept, and set number of identifiers.
      * 
      * \param in  input stream from param file. 
      */
      void readParamType(std::istream& in);

      /**
      * Write type of parameter swept.
      * 
      * \param out  output file stream
      */
      void writeParamType(std::ostream& out) const;

      /**
      * Store the pre-sweep value of the corresponding parameter.
      */
      void getInitial();

      /**
      * Update the corresponding parameter value in the system. 
      * 
      * \param newVal New value of parameter, calculated in specific sweep class.
      */
      void update(double newVal);

      /**
      * Return a string describing the parameter type for this object. 
      */
      std::string type() const;

      /**
      * Accessor for the id_ array.
      * 
      * \param i array index to access
      */
      int id(int i) const
      {  return id_[i];}
      
      /**
      * Return the current system parameter value. 
      */
      double current()
      {  return get_(); }

      /**
      * Return the initial system parameter value. 
      */
      double initial() const
      {  return initial_; }

      /**
      * Return the total change planned for this parameter during sweep.
      */
      double change() const
      {  return change_; }

      /**
      * Set the system associated with this object.
      */
      void setSystem(System<D>& system)
      {  systemPtr_ = &system;}
      
      /**
      * Serialize to or from an archive.
      *
      * \param ar Archive object 
      * \param version archive format version index
      */
      template <class Archive>
      void serialize(Archive ar, const unsigned int version);

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
      
      /// Pointer to the parent system. 
      System<D>* systemPtr_;

      /// Gets the current system parameter value.
      double get_();

      // Set the system parameter value. 
      void set_(double newVal);

   // friends:

      template <int U>
      friend 
      std::istream& operator >> (std::istream&, SweepParameter<U>&);

      template <int U>
      friend 
      std::ostream& 
      operator << (std::ostream&, SweepParameter<U> const&);
   
   };

}
}
// include implementations.
#include "SweepParameter.tpp"

#endif
