#ifndef PRDC_SCFT_REAL_H
#define PRDC_SCFT_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRefReal.h"
#include <util/global.h>
#include <iostream>

namespace Pscf {
namespace Prdc {

   /**
   *
   *
   * \ingroup Pscf_Prdc_Module
   */
   template <int D, class ST>
   class ScftReal : SystemConstRefReal<ST>
   {

   public:

      // Public type name aliases

      /// Base class.
      using Base = SystemConstRefReal<ST>;

      // Inherited type name aliases.
      using SystemT = typename Base::SystemT;
      using MixtureT = typename Base::MixtureT;
      using InteractionT = typename Base::InteractionT;
      using DomainT = typename Base::DomainT;
      using WFieldContainerT = typename Base::WFieldContainerT;
      using CFieldContainerT = typename Base::CFieldContainerT;
      using MaskT = typename Base::MaskT;

      // Public member functions

      /**
      * Constructor.
      */
      ScftReal(SystemT const & system);

      /**
      * Destructor.
      */
      ~ScftReal();

      /**
      * Compute SCFT free energy density and pressure for current fields.
      *
      * Resulting values are retrieved by the fHelmholtz(), fIdeal(), 
      * fInter(), fExt(), and pressure() accessor functions.
      *
      * \pre w().hasData() == true
      * \pre c().hasData() == true
      * \post hasFreeEnergy() == true
      */
      void compute();

      /**
      * Clear all thermodynamic data.
      *
      * \post hasFreeEnergy() == false
      */
      void clear();

      ///@}
      /// \name SCFT Property Accessors
      ///@{

      /**
      * Have free energies and pressure been computed?
      */
      bool hasFreeEnergy() const;

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

      ///@}
      /// \name SCFT Thermodynamic Data Output
      ///@{

      /**
      * Write SCFT thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer, pressure
      * (in units of kT per monomer volume), the volume fraction and
      * chemical potential of each species, and all unit cell parameters.
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

   protected:

      using Base::system;
      using Base::mixture;
      using Base::interaction;
      using Base::domain;
      using Base::c;
      using Base::w;
      using Base::h;
      using Base::mask;
      using Base::fileMaster;

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
      bool hasFreeEnergy_;

   };

   // Inline member functions

   // Get the Helmholtz free energy per monomer / kT.
   template <int D, class ST>
   inline double ScftReal<D,ST>::fHelmholtz() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fHelmholtz_;
   }

   // Get the ideal gas contribution to fHelmholtz.
   template <int D, class ST>
   inline double ScftReal<D,ST>::fIdeal() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fIdeal_;
   }

   // Get the interaction contribution to fHelmholtz.
   template <int D, class ST>
   inline double ScftReal<D,ST>::fInter() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fInter_;
   }

   // Get the external field contribution to fHelmholtz.
   template <int D, class ST>
   inline double ScftReal<D,ST>::fExt() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fExt_;
   }

   // Get the precomputed pressure (units of kT / monomer volume).
   template <int D, class ST>
   inline double ScftReal<D,ST>::pressure() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return pressure_;
   }

   // Have free energies and pressure been computed?
   template <int D, class ST>
   inline bool ScftReal<D,ST>::hasFreeEnergy() const
   {  return hasFreeEnergy_; }

} // namespace Prdc
} // namespace Pscf
#endif
