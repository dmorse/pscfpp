#ifndef PSSP_MIXTURE_TPP
#define PSSP_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <pscf/mesh/Mesh.h>

#include <cmath>

namespace Pscf {
namespace Pssp
{

   template <int D>
   Mixture<D>::Mixture()
    : vMonomer_(1.0),
      ds_(-1.0),
      meshPtr_(0),
      unitCellPtr_(0)
   {  setClassName("Mixture"); }

   template <int D>
   Mixture<D>::~Mixture()
   {}

   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureTmpl< Polymer<D>, Solvent<D> >::readParameters(in);
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);
      read(in, "ds", ds_);

      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
   }

   template <int D>
   void Mixture<D>::setMesh(Mesh<D> const& mesh)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);

      meshPtr_ = &mesh;

      // Set discretization for all blocks
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).setDiscretization(ds_, mesh);
         }
      }

   }

   template <int D>
   void Mixture<D>::setupUnitCell(const UnitCell<D>& unitCell)
   {

      // Set association to unitCell
      unitCellPtr_ = &unitCell;

      for (int i = 0; i < nPolymer(); ++i) {
         polymer(i).setupUnitCell(unitCell);
      }
   }

   /*
   * Compute concentrations (but not total free energy).
   */
   template <int D>
   void Mixture<D>::compute(DArray<Mixture<D>::WField> const & wFields,
                            DArray<Mixture<D>::CField>& cFields)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(mesh().size() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

      int nx = mesh().size();
      int nm = nMonomer();
      int i, j, k;

      // Clear all monomer concentration fields
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == nx);
         UTIL_CHECK(wFields[i].capacity() == nx);
         for (j = 0; j < nx; ++j) {
            cFields[i][j] = 0.0;
         }
      }

      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields);
      }

      // Accumulate monomer concentration fields
      double phi;
      for (i = 0; i < nPolymer(); ++i) {
         phi = polymer(i).phi();
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            int monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            CField& monomerField = cFields[monomerId];
            CField& blockField = polymer(i).block(j).cField();
            for (k = 0; k < nx; ++k) {
               monomerField[k] += phi * blockField[k];
            }
         }
      }

      // To do: Add compute functions and accumulation for solvents.
   }

   /*
   * Compute total stress.
   */
   template <int D>
   void Mixture<D>::computeStress()
   {
      int i, j;

      // Initialize stress to zero
      for (i = 0; i < 6; ++i) {
         stress_[i] = 0;
      }

      // Compute stress for all polymers, after solving MDE
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).computeStress();
      }

      // Accumulate stress for all the polymer chains
      for (i = 0; i < unitCellPtr_->nParameter(); ++i) {
         for (j = 0; j < nPolymer(); ++j) {
            stress_[i] += polymer(j).stress(i);
         }
      }
   }

} // namespace Pssp
} // namespace Pscf
#endif
