#ifndef RP_SCFT_THERMO_H
#define RP_SCFT_THERMO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/SystemConstRef.h>
#include <util/global.h>
#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * %Base class for SCFT thermodynamic property calculators.
   *
   * This class computes and stores values for the SCFT Helmholtz free
   * energy and pressure, and free energy components (ideal, interaction,
   * external). Specializations of this template are used as base classes 
   * for class templates named ScftThermo<D> that are defined in Rpc and 
   * Rpg namespaces.
   *
   * Template parameters:
   *
   *   - D : dimension of space
   *   - T : "Types" class (Rpc::Types<D> or Rpg::Types<D>)
   *
   * \ingroup Rp_Scft_Module
   */
   template <int D, class T>
   class ScftThermo : protected T::SystemConstRef
   {

   public:

      /// \name State modifiers
      ///@{

      /**
      * Compute SCFT free energy density and pressure for current fields.
      *
      * Resulting values are retrieved by the fHelmholtz(), fIdeal(), 
      * fInter(), fExt(), and pressure() accessor functions.
      *
      * \pre w().hasData() == true
      * \pre c().hasData() == true
      * \post hasData() == true
      */
      void compute();

      /**
      * Clear all thermodynamic data.
      *
      * \post hasData() == false
      */
      void clear();

      ///@}
      /// \name SCFT Property Access and Output
      ///@{

      /**
      * Have free energies and pressure been computed?
      */
      bool hasData() const;

      /**
      * Get total Helmholtz free energy per monomer / kT.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fHelmholtz() const;

      /**
      * Get the ideal gas contribution to fHelmholtz.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fIdeal() const;

      /**
      * Get the interaction contribution to fHelmholtz.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fInter() const;

      /**
      * Get the external field contribution to fHelmholtz.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fExt() const;

      /**
      * Get the precomputed pressure times monomer volume / kT.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      * The value is -1 times the grand-canonical free energy per monomer
      * divided by kT.
      */
      double pressure() const;

      /**
      * Write SCFT thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer, pressure
      * (in units of kT per monomer volume), the volume fraction and
      * chemical potential of each species, and all unit cell parameters.
      *
      * If data is not available (i.e,. if hasData() == false), this 
      * function calls compute() on entry.
      *
      * If parameter "out" is a file that already exists, this function
      * will append to the end of that file, rather than overwriting it.
      * Calling writeParamNoSweep and writeThermo in succession with the
      * same output stream will thus produce a single file containing both
      * input parameters and calculated thermodynamic properties.
      *
      * \param out  output stream
      */
      void write(std::ostream& out);

      ///@}

   protected:

      /// Base class type name alias.
      using SystemConstRefT = typename T::SystemConstRef;

      /// Parent System type name alias.
      using SystemT = typename T::System;

      // Protected inherited type name aliases.
      using MixtureT = typename T::Mixture;
      using InteractionT = typename T::Interaction;
      using DomainT = Domain<D,T>;
      using WFieldsT = typename T::WFields;
      using CFieldsT = CFields<D,T>;
      using MaskT = Mask<D,T>;
      using RFieldT = typename T::RField;

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      ScftThermo(SystemT const & system);

      /**
      * Destructor.
      */
      ~ScftThermo();

      // Protected inherited member functions
      using SystemConstRefT::system;
      using SystemConstRefT::mixture;
      using SystemConstRefT::interaction;
      using SystemConstRefT::domain;
      using SystemConstRefT::c;
      using SystemConstRefT::w;
      using SystemConstRefT::h;
      using SystemConstRefT::mask;
      using SystemConstRefT::fileMaster;

   private:

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Ideal gas contribution to fHelmholtz.
      *
      * This includes the internal energy and entropy of
      * non-interacting molecules in the current w fields.
      */
      double fIdeal_;

      /**
      * Interaction contribution to fHelmholtz.
      */
      double fInter_;

      /**
      * External field contribution to fHelmholtz (if any).
      */
      double fExt_;

      /**
      * Pressure times monomer volume / kT.
      *
      * This is -1 times the grand-canonical free energy per monomer,
      * divided by kT.
      */
      double pressure_;

      /**
      * Has SCFT free energy and pressure been computed?
      *
      * This is set true by the compute() member function, and is set 
      * false by the clear() function.
      */
      bool hasData_;

   };

   // Inline member functions

   // Get the Helmholtz free energy per monomer / kT.
   template <int D, class T>
   inline double ScftThermo<D,T>::fHelmholtz() const
   {
      UTIL_CHECK(hasData_);
      return fHelmholtz_;
   }

   // Get the ideal gas contribution to fHelmholtz.
   template <int D, class T>
   inline double ScftThermo<D,T>::fIdeal() const
   {
      UTIL_CHECK(hasData_);
      return fIdeal_;
   }

   // Get the interaction contribution to fHelmholtz.
   template <int D, class T>
   inline double ScftThermo<D,T>::fInter() const
   {
      UTIL_CHECK(hasData_);
      return fInter_;
   }

   // Get the external field contribution to fHelmholtz.
   template <int D, class T>
   inline double ScftThermo<D,T>::fExt() const
   {
      UTIL_CHECK(hasData_);
      return fExt_;
   }

   // Get the precomputed pressure (units of kT / monomer volume).
   template <int D, class T>
   inline double ScftThermo<D,T>::pressure() const
   {
      UTIL_CHECK(hasData_);
      return pressure_;
   }

   // Have free energies and pressure been computed?
   template <int D, class T>
   inline bool ScftThermo<D,T>::hasData() const
   {  return hasData_; }

} // namespace Rp
} // namespace Pscf
#endif
