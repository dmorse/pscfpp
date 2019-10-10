#ifndef PSPG_POLYMER_H
#define PSPG_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/solvers/PolymerTmpl.h>
#include <pspg/field/RDField.h>
#include <util/containers/FArray.h> 

namespace Pscf { 
namespace Pspg { 


   /**
   * Descriptor and solver for a branched polymer species.
   *
   * The block concentrations stored in the constituent Block<D>
   * objects contain the block concentrations (i.e., volume 
   * fractions) computed in the most recent call of the compute 
   * function.
   *
   * The phi() and mu() accessor functions, which are inherited 
   * from PolymerTmp< Block<D> >, return the value of phi (spatial 
   * average volume fraction of a species) or mu (chemical
   * potential) computed in the last call of the compute function.
   * If the ensemble for this species is closed, phi is read from 
   * the parameter file and mu is computed. If the ensemble is
   * open, mu is read from the parameter file and phi is computed.
   *
   * \ingroup Pspg_Solvers_Module
   */
   template <int D>
   class Polymer : public PolymerTmpl< Block<D> >
   {

   public:

      typedef PolymerTmpl< Block<D> > Base;

      typedef typename Block<D>::WField  WField;

      Polymer();

      ~Polymer();

      void setPhi(double phi);

      void setMu(double mu);

      /**
      * Compute solution to MDE and concentrations.
      */ 
      void setupUnitCell(UnitCell<D> const & unitCell);

      /**
      * Compute solution to MDE and concentrations.
      */ 
      void compute(DArray<WField> const & wFields);

      // Stress due to a whole Polymer chain
      FArray<double, 6> PcStress;
      
      /**
      * Compute stress from a polymer chain, needs a pointer to basis
      */
      void ComputePcStress(WaveList<D>& wavelist);

      using Base::nBlock;
      using Base::block;
      using Base::ensemble;
      using Base::solve;
      using Base::length;

   protected:

      using ParamComposite::setClassName;

   private: 

      int nParams_;
      using Base::phi_;
      using Base::mu_;

   };

}
}
#include "Polymer.tpp"
#endif
