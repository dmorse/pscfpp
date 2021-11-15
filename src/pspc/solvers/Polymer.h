#ifndef PSPC_POLYMER_H
#define PSPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/solvers/PolymerTmpl.h>
#include <pspc/field/RField.h>
#include <util/containers/FArray.h>      // member template

namespace Pscf { 
namespace Pspc { 

   /**
   * Descriptor and solver for one polymer species.
   *
   * The phi() and mu() accessor functions, which are inherited 
   * from PolymerTmp< Block<D> >, return the value of phi (spatial 
   * average volume fraction of a species) or mu (chemical
   * potential) computed in the last call of the compute function.
   * If the ensemble for this species is closed, phi is read from 
   * the parameter file and mu is computed. If the ensemble is
   * open, mu is read from the parameter file and phi is computed.
   *
   * The block concentrations stored in the constituent Block<D> objects
   * contain the block concentrations (i.e., volume fractions) computed 
   * in the most recent call of the compute function. These can be 
   * accessed using the Block<D>::cField()
   *
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Polymer : public PolymerTmpl< Block<D> >
   {

   public:

      /**
      * Base class typedef (PolymerTmpl instance)
      */
      typedef PolymerTmpl< Block<D> > Base;

      /**
      * Chemical potential field typedef
      */
      typedef typename Block<D>::WField  WField;

      /**
      * Default constructor.
      */
      Polymer();

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Set value of phi (volume fraction), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param phi  desired volume fraction for this species
      */
      void setPhi(double phi);

      /**
      * Set value of mu (chemical potential), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param mu  desired chemical potential for this species
      */
      void setMu(double mu);

      /**
      * Set up the unit cell after a change in unit cell parameters.
      *
      * This function should be called after each change in the unit
      * cell. It sets up unit cell information for all blocks in this
      * polymer.
      *
      * \param unitCell crystallographic unit cell
      */ 
      void setupUnitCell(UnitCell<D> const & unitCell);

      /**
      * Compute solution to MDE and block concentrations.
      * 
      * This function sets up w-fields in the MDE solvers for all blocks
      * and then calls the base class PolymerTmpl solve function. This
      * base class solve function solves the MDE for all propagators
      * propagators and computes the properly scaled volume fraction 
      * fields for all blocks. After this function is called, the
      * associated Block objects store pre-computed propagator 
      * solutions and block volume fraction fields. 
      *
      * \param wFields array of chemical potential fields.
      */ 
      void compute(DArray<WField> const & wFields);

      /**
      * Compute stress from a polymer chain.
      *
      * This function computes contributions from this species to the 
      * derivatives of free energy per monomer with respect to unit cell 
      * parameters and stores the values.
      */
      void computeStress();

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *  
      * This function gets the precomputed value of the derivative of
      * free energy per monomer with respect to unit cell parameter n.
      * The value must have been pre-computed by calling computeStress().
      *  
      * \param n index of unit cell parameter
      */
      double stress(int n) const;

      // Inherited public functions

      using Base::nBlock;
      using Base::block;
      using Base::ensemble;
      using Base::solve;
      using Base::length;

   protected:

      using ParamComposite::setClassName;
      using Base::phi_;
      using Base::mu_;

   private: 

      /// Stress contribution from this polymer species
      FArray<double, 6> stress_;

      /// Pointer to associated UnitCell<D>
      const UnitCell<D>* unitCellPtr_;

   };

   /// Stress with respect to unit cell parameter n.
   template <int D>
   inline double Polymer<D>::stress(int n) const
   {  return stress_[n]; }
  
   #ifndef PSPC_POLYMER_TPP
   // Supress implicit instantiation
   extern template class Polymer<1>;
   extern template class Polymer<2>;
   extern template class Polymer<3>;
   #endif

}
}
#endif
