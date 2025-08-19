#ifndef PRDC_SCFT_THERMO_TMPL_H
#define PRDC_SCFT_THERMO_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRefTmpl.h"
#include <util/global.h>
#include <iostream>

namespace Pscf {
namespace Prdc {

   /**
   * %Base class for SCFT thermodynamic property calculators.
   *
   * This class computes and stores values for the SCFT Helmholtz
   * free energy and pressure, and free energy components (ideal,
   * interaction, external). It is used as a base class for 
   * classes named ScftThermo defined in Rpc and Rpg namespaces.
   *
   * \ingroup Prdc_System_Module
   */
   template <int D, class ST>
   class ScftThermoTmpl : protected SystemConstRefTmpl<ST>
   {

   public:

      /// Base class type name alias.
      using Base = SystemConstRefTmpl<ST>;

      /// Parent System type name alias.
      using SystemT = typename Base::SystemT;

      // Public member functions

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      ScftThermoTmpl(SystemT const & system);

      /**
      * Destructor.
      */
      virtual ~ScftThermoTmpl();

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
      * If data is not available (i.e,. if hasData() == false), this function 
      * calls compute() on entry.
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

      // Protected inherited type name aliases.
      using MixtureT = typename Base::MixtureT;
      using InteractionT = typename Base::InteractionT;
      using DomainT = typename Base::DomainT;
      using WFieldsT = typename Base::WFieldsT;
      using CFieldContainerT = typename Base::CFieldContainerT;
      using MaskT = typename Base::MaskT;
      using RFieldT = typename Base::RFieldT;

      // Protected inherited member functions
      using Base::system;
      using Base::mixture;
      using Base::interaction;
      using Base::domain;
      using Base::c;
      using Base::w;
      using Base::h;
      using Base::mask;
      using Base::fileMaster;

      /**
      * Inner product of two fields.
      *
      * \param A  1st field
      * \param B  2nd field
      * \return inner product (sum of product A[i]*B[i] over mesh points)
      */
      virtual
      double innerProduct(RFieldT const & A, RFieldT const & B) const
      {  return 0.0; };

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
      * This is set true in the compute() function, and is set false by
      * clear function.
      */
      bool hasData_;

   };

   // Inline member functions

   // Get the Helmholtz free energy per monomer / kT.
   template <int D, class ST>
   inline double ScftThermoTmpl<D,ST>::fHelmholtz() const
   {
      UTIL_CHECK(hasData_);
      return fHelmholtz_;
   }

   // Get the ideal gas contribution to fHelmholtz.
   template <int D, class ST>
   inline double ScftThermoTmpl<D,ST>::fIdeal() const
   {
      UTIL_CHECK(hasData_);
      return fIdeal_;
   }

   // Get the interaction contribution to fHelmholtz.
   template <int D, class ST>
   inline double ScftThermoTmpl<D,ST>::fInter() const
   {
      UTIL_CHECK(hasData_);
      return fInter_;
   }

   // Get the external field contribution to fHelmholtz.
   template <int D, class ST>
   inline double ScftThermoTmpl<D,ST>::fExt() const
   {
      UTIL_CHECK(hasData_);
      return fExt_;
   }

   // Get the precomputed pressure (units of kT / monomer volume).
   template <int D, class ST>
   inline double ScftThermoTmpl<D,ST>::pressure() const
   {
      UTIL_CHECK(hasData_);
      return pressure_;
   }

   // Have free energies and pressure been computed?
   template <int D, class ST>
   inline bool ScftThermoTmpl<D,ST>::hasData() const
   {  return hasData_; }

} // namespace Prdc
} // namespace Pscf
#endif
