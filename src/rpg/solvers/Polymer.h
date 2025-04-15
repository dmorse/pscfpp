#ifndef RPG_POLYMER_H
#define RPG_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PolymerTmpl.h>      // base class template
#include <rpg/solvers/Block.h>             // base class argument

//#include <prdc/cuda/RField.h>
#include <util/containers/FArray.h> 

// Forward declarations
namespace Util {
   template <typename T> class DArray;
}
namespace Pscf { 
   namespace Prdc {
      template <int D> class UnitCell;
      namespace Cuda {
         template <int D> class RField;
      }
   }
   namespace Rpg { 
      template <int D> class WaveList;
   }
}

namespace Pscf { 
namespace Rpg { 

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Descriptor and solver for a branched polymer species.
   *
   * The propagator solutions and block concentrations stored in the
   * constituent Block<D> objects contain solutions computed in the most
   * recent call of the compute function.
   *
   * The phi() and mu() accessor functions, which are inherited from
   * PolymerSpecies, return the value of phi (spatial average volume
   * fraction of a species) or mu (species chemical potential) computed 
   * in the last call of the compute function.  If the ensemble for this 
   * species is closed, phi is read from the parameter file and mu is 
   * computed. If the ensemble is open, mu is read from the parameter 
   * file and phi is computed.
   *
   * \ingroup Rpg_Solvers_Module
   */
   template <int D>
   class Polymer : public PolymerTmpl< Block<D> >
   {

   public:

      /**
      * Alias for base class.
      */
      typedef PolymerTmpl< Block<D> > Base;

      /**
      * Constructor. 
      */
      Polymer();

      /**
      * Destructor. 
      */
      ~Polymer();

      /**
      * Set the phi (volume fraction) for this species.
      *
      * \param phi volume fraction (input)
      */
      void setPhi(double phi);

      /**
      * Set the mu (chemical potential) for this species.
      * 
      * \param mu chemical potential (input)
      */
      void setMu(double mu);

      /**
      * Store the number of lattice parameters in the unit cell.
      *
      * \param nParams  the number of lattice parameters in the unit cell
      */ 
      void setNParams(int nParams);

      /**
      * Compute MDE solutions and block concentrations.
      * 
      * This function sets up w-fields in the MDE solvers for all blocks
      * and then calls the base class PolymerTmpl solve function. This
      * solves the MDE for all propagators and computes the properly 
      * scaled volume fraction fields for all blocks. After this function 
      * is called, the associated Block objects store pre-computed 
      * propagator solutions and block volume fraction fields. 
      *
      * The parameter phiTot is only relevant to problems such as thin
      * films in which the material is excluded from part of the unit
      * cell by imposing an inhogeneous constraint on the sum of the
      * monomer concentrations (i.e., a "mask"). 
      *
      * \param wFields array of chemical potential fields.
      * \param phiTot  volume fraction of unit cell occupied by material
      */ 
      void compute(DArray< RField<D> > const & wFields, 
                   double phiTot = 1.0);

      /**
      * Compute stress contribution from this polymer species.
      */
      void computeStress();

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * Get the contribution from this polymer species to the derivative of
      * free energy per monomer with respect to unit cell parameter n.
      *
      * \param n unit cell parameter index
      */
      double stress(int n);

      // Inherited public functions
      using Base::edge;
      using Base::block;
      using Base::propagator;
      using PolymerSpecies::vertex;
      using PolymerSpecies::propagatorId;
      using PolymerSpecies::path;
      using PolymerSpecies::nBlock;
      using PolymerSpecies::nVertex;
      using PolymerSpecies::nPropagator;
      using PolymerSpecies::length;
      using PolymerSpecies::nBead;
      using PolymerSpecies::type;
      using Species::phi;
      using Species::mu;
      using Species::q;
      using Species::ensemble;

   protected:

      // Protected inherited function with non-dependent names
      using ParamComposite::setClassName;

   private: 

      /// Stress contribution from this polymer species.
      FArray<double, 6> stress_;
     
      /// Number of unit cell parameters. 
      int nParams_;

      // Restrict access
      using Base::solve;
      using Species::phi_;
      using Species::mu_;
      using Species::q_;
      using Species::ensemble_;

   };

   template <int D>
   double Polymer<D>::stress(int n)
   {  return stress_[n]; }

   #ifndef RPG_POLYMER_TPP
   // Suppress implicit instantiation
   extern template class Polymer<1>;
   extern template class Polymer<2>;
   extern template class Polymer<3>;
   #endif

}
}
//#include "Polymer.tpp"
#endif


