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
	template <int D>
   Basis<D>::Basis()
   	: nWave_(0), nStar_(0), unitCellPtr_(0), mesh_(0)
   {}

   template <int D>
   void Basis<D>::makeBasis(const Mesh<D>& mesh, const UnitCell<D>& unitCell,
                            std::string groupName)
   {
   	if(groupName == "I") {
         //To Do: Need some method to obtain the space group symmetry 
         //       if other than I
   	}
   	else{
   		UTIL_THROW("Unimplemented space group");
   	}

      //Associate both mesh and unit cell
      mesh_ = &mesh;
      unitCellPtr_ = &unitCell;
      
      meshIterator<D> itr(mesh.dimensions());
      
      //make waves
      waves_.allocate(mesh.size());
      waveId_.allocate(mesh.size());
      nWave_ = mesh.size();

      //requires the right implementation of itr.position(). Fixed in 
      //requests mde pull
      for(itr.begin(); !itr.atEnd(); ++itr) {
         waves_[itr.rank()].sqNorm = unitCell.ksq(itr.position());
         waves_[itr.rank()].indicesDft = itr.position();

         if(itr.position(D-1) > mesh_->dimension(D-1)/2 + 1) {
            waves_[itr.rank()].implicit = true;
         }

         //`unsorted' waves; waves appear in grid order
         waveId_[itr.rank()] = itr.rank();
      }
      
      //make stars
      //To do: If not I, sort according to increasing ksq values
      //Changed: The size of stars_ is not known apriori. Used a GArray
      //For `I', reserved the right size.
      stars_.reserve(mesh.size());
      int beginWave = 0;
      int endWave = 0;
      int iStar = 0;
      bool cancel;
      int invertFlag = 0;
      std::complex<double> cNorm;

      for (int i = 1; i < mesh.size(); i++) {
         //To do: If not I, create a list of waves with equal ksq 
         //magnitudes. The I space group implicitly assumes each waves belong
         //in a single star and the sorting is not done.
         
         endWave = beginWave+1;
         double starPhase[1] = {0};
         //To do: unimplemented local scope function
         //cancel = isCancelled(waves_[beginWave].indices);
         cancel = false;

         #if 0
         //check for cancellation consistency in star
         for(int j = beginWave; j < endWave; ++j) {
            if( cancel != isCancelled(waves_[j].indices))
               UTIL_THROW("Inconsistent cancellation in star");
         }
         #endif

         if(!cancel) {
            if(invertFlag == 1)
            { 
               cNorm = exp(std::complex<double>(0,1) *
                          starPhase[endWave-1-beginWave]);
            }
            else 
            { cNorm = exp(std::complex<double>(0,1)*starPhase[0]);}

            cNorm *= sqrt(double(endWave - beginWave));
         }

         for(int j = beginWave; j < endWave; ++j) {
            waves_[j].starId = iStar;
            if(!cancel)
            {
               waves_[j].coeff = exp(std::complex<double>(0,1) *
                                    starPhase[j - beginWave])/cNorm;
            }
            else
            {   waves_[j].coeff = exp(std::complex<double>(0,0)); }
         }
         
         //fill up stars_ object
         //unimplemented: sign flag.
         stars_[iStar].size = endWave - beginWave;
         stars_[iStar].beginId = beginWave;
         stars_[iStar].endId = endWave - 1;
         stars_[iStar].cancel = cancel;

         //To do: a method here to decide what is the actual invertFlag
         //A lot easier when waves are sorted in min(ksq). Currently, the second
         //wave of 001 index 1 is paired with say 006. So we default to 0
         stars_[iStar].invertFlag = invertFlag;

         if(invertFlag == -1)
            stars_[iStar].waveBz = shiftToMinimum(waves_[endWave-1].indicesDft,
                                      mesh_->dimensions(), unitCellPtr_);
         else
            stars_[iStar].waveBz = shiftToMinimum(waves_[beginWave].indicesDft,
                                      mesh_->dimensions(), unitCellPtr_);

         iStar++;
         beginWave = endWave;
      }
   }


   template <int D>
   void Basis<D>::convertFieldComponentsToDFT(DArray<double>& components, RFieldDft<D>& dft)
   {
      std::complex<double> coeff;
      double z;
      
      //loop through a factor of 2 more indices. 
      //Difficult to optimize without specialization or making new objects
      for(int i = 0; i < mesh_->size(); i++){
         
         //if last indice of waves is not > n3/2+1
         int waveId = wave(mesh_->position[i]);
         int starId = waves_[waveId].starId;
         if(!waves_[waveId].implicit)
         {
            coeff = waves_[waveId].coeff;
            //determining the rank of RFieldDft
            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if(j == D-1){
                  offsets[j-1] = offsets[j]*(mesh_->dimension(j)/2 + 1);
               }
               offsets[j-1] = offsets[j]*mesh_->dimensions(j);
            }
            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += waves_[waveId].indicesDft[j] * offsets[j];
            }


            switch(stars_[starId].invertFlag){
               case(0) :
                  z = components[starId];
                  dft[rank] = z * coeff;
               case(1) :
                  UTIL_THROW("Not an implemented function");
                  //z1 z2 not implemented correctly
                  //z1 = components[starId];
                  //z2 = components[starId+1];
                  //dft[rank] = fftw_complex(z1, -z2) * coeff / sqrt(2);
               case(-1) :
                  UTIL_THROW("Not an implemented function");
                  //z1 z2 not implemented correctly
                  //z1 = components[starId-1];
                  //z2 = components[starId];
                  //dft[rank] = fftw_complex(z1, z2) * coeff /sqrt(2);
            }

         }
         else{
            //do nothing. Dft does not need allocation for this
         }
      }



   }

   template <int D>
   void Basis<D>::convertFieldDftToComponents(RFieldDft<D>& dft, DArray<double>& components)
   {
      //if its not too distrupting maybe use consistent names for logical size?
      int nStars = stars_.size();
      UTIL_CHECK(nStars == components.capacity());
      IntVec<D> indiceBz;
      std::complex<double> coeff;
      std::complex<double> temp;
      fftw_complex z;
      for(int i = 0; i < nStars; ++i) {
         indiceBz = stars_[i].waveBz;
         coeff = waves_[wave(indiceBz)].coeff;

         //check if wave exists in dft
         if(!waves_[wave(indiceBz)].implicit) {
            mesh_->shift(indiceBz);
            //might be good to rethink how elements of RFieldDft is accessed
            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if(j == D-1){
                  offsets[j-1] = offsets[j]*(mesh_->dimension(j)/2 + 1);
               }
               offsets[j-1] = offsets[j]*mesh_->dimensions(j);
            }
            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }
            z = dft[rank];
            //Note: the code implies that we do not need indicesDft from Wave? 
         }
         else{
         //wave does not exists in dft. have to find explicit pair
            mesh_->shift(indiceBz);
            indiceBz[D-1] = -indiceBz[D-1];
            mesh_->shift(indiceBz);

            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if(j == D-1){
                  offsets[j-1] = offsets[j]*(mesh_->dimension(j)/2 + 1);
               }
               offsets[j-1] = offsets[j]*mesh_->dimensions(j);
            }
            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }
            z = dft[rank];
            z[1] = -z[1];
         }
         
         //reintepret cast is not used since coding standards is old
         temp = std::complex<double>(z[0],z[1]);
         //assign value to components
         switch(stars_[i].invertFlag) {
            case 0 :
               components[i] = (temp/coeff).real();
            case 1 :
               UTIL_THROW("Not an implemented function");
               //components[i] = double(temp/coeff) * sqrt(2);
            case -1 :
               UTIL_THROW("Not an implemented function");
               //components[i] = imag(temp/coeff) * sqrt(2);
         }
      }
   }

   template <int D>
   int Basis<D>::nBasis() const
   {
      int count = 0;
      for (int i = 0; i < stars_.capacity(); ++i) {
         if (!stars_[i].cancel)
            count++;
      }
      return count;
   }

   template <int D>
   typename Basis<D>::Wave& Basis<D>::wave(IntVec<D> vector)
   {
      if(!mesh_->isInMesh(vector))
         mesh_->shift(vector);
   	return waves_[waveId_[mesh_->rank(vector)]];
   }


}
}

#endif