#ifndef PSSP_BASIS_TPP
#define PSSP_BASIS_TPP
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Basis.h"

namespace Pscf {
namespace Pssp
{
   using namespace Util;

   template <int D>
   Basis<D>::Basis()
    : nWave_(0),
      nStar_(0),
   {
   	//struct cannot have constructors
   	unitCell_ = 0;
   	mesh_ = 0;
   }

   template <int D>
   void Basis<D>::setUnitCell(const UnitCell<D>& unitCell)
   { unitCell = &(unitCell); }

   template <int D>
   void Basis<D>::setMesh(const Mesh<D>& mesh)
   { mesh_ = &(mesh); }

   template <int D>
   void Basis<D>::allocate()
   {
   	UTIL_CHECK( mesh_ != 0);
   	UTIL_CHECK( unitCell_ != 0);

   	//only true for space group 1
   	//check for cancellation otherwise
   	int size = mesh_->size();
   	waves_.allocate(size);
   	nWave_ = size;

   	//only true for space group 1
   	//check for cancellation otherwise
   	stars_.allocate(size);
   	nStar_ = size;
   }

   template <int D>
   void Basis<D>::deallocate()
   {
   	waves_.deallocate();
   	stars_.deallocate();
   }

   template <int D>
   int Basis<D>::nWave() const
   { return nWave_; }

   template <int D>
   int Basis<D>::nStar() const
   { return nStar_; }

   template <int D>
   Basis<D>::Wave& Basis<D>::wave(int i)
   { return waves_[i]; }

   template <int D>
   Basis<D>::Wave& wave(IntVec<D> vector)
   { 
   	UTIL_CHECK(mesh_->isInMesh(vector));
   	/*int rank;
   	for (int i = 0; i < D; i++){
         i = vector[0] * vector[1] * vector[2];
   	}
   	return waves_[i]; */
   }
}
}

#endif