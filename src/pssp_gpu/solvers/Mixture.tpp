#ifndef PSSP_GPU_MIXTURE_TPP
#define PSSP_GPU_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <pssp_gpu/GpuResources.h>
//#include <fd1d/domain/Domain.h>

#include <cmath>

//in propagator.h
//static __global__ void assignUniformReal(cufftReal* result, cufftReal uniform, int size);

//theres a precision mismatch here. need to cast properly.
static __global__ void accumulateConc(cufftReal* result, double uniform, cufftReal* cField, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] += (float)uniform * cField[i];
   }
}

namespace Pscf { 
namespace Pssp_gpu
{ 

   template <int D>
   Mixture<D>::Mixture()
    : vMonomer_(1.0),
      ds_(-1.0),
      meshPtr_(0)
   {  setClassName("Mixture"); }

   template <int D>
   Mixture<D>::~Mixture()
   {}

   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureTmpl< Pscf::Pssp_gpu::Polymer<D>, Pscf::Pssp_gpu::Solvent<D> >::readParameters(in);
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
      nParams_ = unitCell.nParams();
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
      int i, j;

      // Clear all monomer concentration fields
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == nx);
         UTIL_CHECK(wFields[i].capacity() == nx);
         //cFields[i][j] = 0.0;
         assignUniformReal<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(cFields[i].cDField(), 0.0, nx);
      }

      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields);
      }

      // Accumulate monomer concentration fields
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            int monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            CField& monomerField = cFields[monomerId];
            CField& blockField = polymer(i).block(j).cField();
            //monomerField[k] += polymer(i).phi() * blockField[k];
            accumulateConc<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(monomerField.cDField(), 
                        polymer(i).phi(), blockField.cDField(), nx);
         }
      }
      // To do: Add compute functions and accumulation for solvents.

   }

   /*  
   * Compute Total Stress.
   */  
   template <int D>
   void Mixture<D>::computeTStress(WaveList<D>& wavelist)
   {   
      int i, j;

      for (i = 0; i < nParams_; ++i) {
         TStress [i] = 0;
         //std::cout<<"TStress ["<<i<<"] = "<< TStress [i]<<std::endl;
      } 
      // Compute Stress for all polymers, must be done after computing concentrations

      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).ComputePcStress(wavelist);
      } 
      std::cout<<"---------"<<std::endl;
      // Accumulate stress for all the polymer chains
      //why is done over all 6 parameter?
      for (i = 0; i < nParams_; ++i) {
         for (j = 0; j < nPolymer(); ++j) {
            TStress [i] += polymer(j).PcStress[i];
          //  std::cout<<"PcStress ["<<j<<","<<i<<" ] = "<< polymer(j).PcStress[i] <<std::endl;
           // std::cout<<"TStress ["<<i<<"] = "<< TStress [i]<<std::endl;
         }   
      }   
   }  

} // namespace Pssp_gpu
} // namespace Pscf


#endif
